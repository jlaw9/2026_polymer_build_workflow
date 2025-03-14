'''
The structure of a polymer building project workflow,
including a shared namespace, labels, conditions, and operations
'''

import logging
import warnings
warnings.catch_warnings(record=True)

import sys, time
from argparse import ArgumentParser, Namespace

from typing import ClassVar, Iterable, Optional, Union

import pickle, json
from functools import partial
from pathlib import Path

import numpy as np

# OpenFF toolkits
from openff.toolkit import Molecule, Topology, ForceField
from openff.toolkit.utils.exceptions import (
    MoleculeParseError,
    UnassignedChemistryInPDBError,
    IncorrectNumConformersWarning,
    InconsistentStereochemistryError,
    ToolkitUnavailableException,
)
logging.getLogger('openff.toolkit.typing.engines.smirnoff.parameters').setLevel(logging.CRITICAL) # silence annoying Electrostatics up-conversion INFO logs
warnings.filterwarnings('ignore', category=IncorrectNumConformersWarning) # silence annoying Conformers warning from Espaloma

from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper
try: # attempt to deregister OpenEye toolkit
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper) # avoid expensive and non-standard OpenEye operations
except ToolkitUnavailableException:
    pass

from openff.units import (
    unit as offunit,
    Unit as OFFUnit,
    Quantity as OFFQuantity,
)

from openff.interchange import Interchange

# Cheminformatics
from rdkit import Chem
from rdkit.Chem.rdmolops import AromaticityModel, SanitizeFlags
from rdkit.Chem.rdqueries import XAtomQueryAtom, MAtomQueryAtom, AtomNumEqualsQueryAtom

# OpenMM imports
from openmm import LangevinMiddleIntegrator
from openmm import Context, XmlSerializer

from openmm.unit import (
    Quantity as OMMQuantity,
    Unit as OMMUnit,
)
from openmm.unit import femtosecond, picosecond, kelvin, kilojoule_per_mole

# Signac imports
from signac.job import Job
from flow import FlowProject

# Custom (polymerist) imports
import polymerist as ps
from polymerist.genutils.logutils.IOHandlers import submodule_loggers, get_active_loggers
POLYMERIST_LOGGERS = [
    logger
        for logger in submodule_loggers(ps).values()
            if (logger is not None) and (not isinstance(logger, logging.PlaceHolder))
] # TODO: move this into polymerist?

from polymerist.polymers.monomers import MonomerGroup
from polymerist.polymers.building import build_linear_polymer, mbmol_to_openmm_pdb

from polymerist.mdtools.openfftools import topology, boxvectors
from polymerist.mdtools.openfftools.partition import partition
from polymerist.mdtools.openfftools.unitsys import openmm_to_openff
from polymerist.mdtools.openfftools.partialcharge.molchargers import MolCharger

from polymerist.mdtools.openmmtools.serialization import apply_state_to_context
from polymerist.mdtools.openmmtools.evaluation import eval_openmm_energies_separated
from polymerist.mdtools.lammpstools.lammpseval import get_lammps_energies

from polymerist.rdutils.rdcoords.tiling import rdmol_effective_radius
from polymerist.rdutils.reactions.reactions import AnnotatedReaction, BadNumberReactants
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.substructures import num_substruct_queries_distinct
from polymerist.rdutils.sanitization import sanitize_mol

from polymerist.smileslib.cleanup import expanded_SMILES

# Utils imports - made these non-relative to avoid screwing up external vs internal call
try: # call as python module
    from .utils.logs import redirect_to_logfile
    from .utils.filelib import is_empty
    from .utils.cheminf import generate_smarts_fragments
    from .utils.dataIO import read_rxn_mapping_data
    from .utils.offlib import elem_counts
    from .utils.packing import generate_uniform_subpopulated_lattice
    from .utils.mdexport import interchange_to_openmm
    from .reactions import RXNS_DIR
    from .environments.cuboulder import CUAlpineEnvironment, CUBlancaShirtsEnvironment # inject CURC-specific environment config
except ImportError: # call as script file
    from utils.logs import redirect_to_logfile
    from utils.filelib import is_empty
    from utils.cheminf import generate_smarts_fragments
    from utils.dataIO import read_rxn_mapping_data
    from utils.offlib import elem_counts
    from utils.packing import generate_uniform_subpopulated_lattice
    from utils.mdexport import interchange_to_openmm
    from reactions import RXNS_DIR
    from environments.cuboulder import CUAlpineEnvironment, CUBlancaShirtsEnvironment # inject CURC-specific environment config


# DEFINING THE SIGNAC PROJECT CLASS PROPER 
class PolymerBuildProject(FlowProject):
    '''Project for automated high-throughput generation of polymer structure and MD inputs from chemical data'''
    # GLOBAL CONFIG - TODO: make these configuratble via argparse to the containing script
    ## LOGGING AND REPORTING FORMATS
    QUANTITY_PRECISION  : ClassVar[int] = 4 # number of decimal place to report Quantities when logging
    LOGLEVEL            : ClassVar[int] = logging.INFO
    ENERGY_UNIT         : ClassVar[OMMUnit] = kilojoule_per_mole 
    OP_TIME_RECORD_NAME : ClassVar[str] = 'operation_times_sec'
    
    ## CHEMICAL SPECIFICATION
    RELAXED_STEREO      : ClassVar[bool] = True
    N_ATOM_CAP_MONOMER  : ClassVar[int] = 150
    SANITIZE_OPS        : ClassVar[SanitizeFlags] = SanitizeFlags.SANITIZE_ALL
    AROMATICITY_MODEL   : ClassVar[AromaticityModel] = AromaticityModel.AROMATICITY_MDL
    REGISTERED_RXN_SMARTS : ClassVar[dict[str, str]] = {}
    
    BLACKLISTED_ATOM_QUERIES : ClassVar[dict[str, Chem.QueryAtom]] = {
        'boron'      : AtomNumEqualsQueryAtom(5, negate=False),
        'silicon'    : AtomNumEqualsQueryAtom(14, negate=False),
        # 'phosphorus' : AtomNumEqualsQueryAtom(15, negate=False),
        'sulfur'     : AtomNumEqualsQueryAtom(16, negate=False),
        'metal'      : MAtomQueryAtom(),
        # 'halogen'    : XAtomQueryAtom(),
    }

    # PROJECT-WIDE FILE NAMES
    ## STRUCTURE FILES
    LOGFILE_NAME        : ClassVar[str] = 'build_logs.log'
    FRAGMENTS_PATH      : ClassVar[str] = 'fragments.json'
    OLIGOMER_PDB        : ClassVar[str] = 'oligomer.pdb'
    OLIGOMER_SDF        : ClassVar[str] = 'oligomer.sdf'
    MELT_NEAT_SDF       : ClassVar[str] = 'melt_neat.sdf'
    INTERCHANGE_PATH    : ClassVar[str] = 'interchange.pkl'

    # LAMMPS
    LAMMPS_DIR          : ClassVar[str] = 'LAMMPS'
    LAMMPS_INPUT_PATH   : ClassVar[str] = f'{LAMMPS_DIR}/inputs.in'
    LAMMPS_DATA_PATH    : ClassVar[str] = f'{LAMMPS_DIR}/data.lmp'
    LAMMPS_PATHS : ClassVar[list[str]] = (
        LAMMPS_INPUT_PATH,
        LAMMPS_DATA_PATH,
    )
    LAMMPS_ENERGIES     : ClassVar[str] = f'{LAMMPS_DIR}/energies_lammps.json' # NOTE: deliberately NOT be lumped w/ MD input files

    ## OpenMM
    OPENMM_DIR          : ClassVar[str] = 'OpenMM'
    OPENMM_STATE_PATH   : ClassVar[str] = f'{OPENMM_DIR}/state.xml'
    OPENMM_SYSTEM_PATH  : ClassVar[str] = f'{OPENMM_DIR}/system.xml'
    OPENMM_TOPO_PATH    : ClassVar[str] = f'{OPENMM_DIR}/topology.pdb'
    OPENMM_INTEG_PATH   : ClassVar[str] = f'{OPENMM_DIR}/integrator.xml'
    OPENMM_PATHS : ClassVar[list[str]] = (
        OPENMM_STATE_PATH,
        OPENMM_SYSTEM_PATH,
        OPENMM_TOPO_PATH,
        OPENMM_INTEG_PATH,
    )
    OPENMM_ENERGIES     : ClassVar[str] = f'{OPENMM_DIR}/energies_openmm.json' # NOTE: deliberately NOT be lumped w/ MD input files
    
# CHEMICAL VALIDATION AND PERCEPTION UTILITIES - !!MOVE THESE TO UTILS SUBMODULE EVENTUALLY!!
def sanitized_mol_from_smiles(
        smiles : str,
        separate_mols : bool=True,
        sanitize_ops : SanitizeFlags=PolymerBuildProject.SANITIZE_OPS,
        aromaticity_model : AromaticityModel=PolymerBuildProject.AROMATICITY_MODEL,
    ) -> Union[Chem.Mol, tuple[Chem.Mol]]:
    '''Load a mol from SMILES and apply the predefined sanitization and aromaticity operations'''
    mol = Chem.MolFromSmiles(smiles, sanitize=False) # CRITICAL that sanitize=False to avoid stripping
    sanitize_mol(mol, sanitize_ops=sanitize_ops, aromaticity_model=aromaticity_model, in_place=True)
    
    if separate_mols:
        return Chem.GetMolFrags(mol, asMols=True)
    return mol

## ATOMS, MONOMERS, AND REACTION MECHANISMS WHICH ARE, FOR ONE REASON OR ANOTHER, NOT ALLOWED
# DEVNOTE: one would think this could be automated by setting a cap on the number of automorphisms, but this cap grows far too quickly
# with the computational cost of evaluating those automorhpisms (via cap on number of substruct matches) to be work it
_blacklisted_monomer_smiles = [ # monomers which are, for one reason or another, disallowed
    'CC(C)(C)c1cc(c(Oc2ccc(cc2)N(c3ccc(N)cc3)c4ccc(N)cc4)c(c1)C(C)(C)C)C(C)(C)C',  # the extraordinary number of symmetries of this amine ("4-N-(4-aminophenyl)-4-N-[4-(2,4,6-tritert-butylphenoxy)phenyl]benzene-1,4-diamine")... 
    'CC(C)(C)c1cc(Oc2ccc(-c3ccc(N)cc3)cc2C(F)(F)F)c(C(C)(C)C)cc1Oc1ccc(-c2ccc(N)cc2)cc1C(F)(F)F', # ...mean it takes impractically long to isomorphism match during the Topology partition step
    'CCCCCCCCCCCCCCCCC(CO)C(CO)CCCCCCCCCCCCCCCC',
] 
BLACKLISTED_MONOMER_MOLS : dict[str, Chem.Mol] = {}
for smiles in _blacklisted_monomer_smiles:
    exp_smi = expanded_SMILES(smiles, assign_map_nums=False)
    BLACKLISTED_MONOMER_MOLS[smiles] = sanitized_mol_from_smiles(exp_smi, separate_mols=False) # though single-molecules, need to separate to avoid tuple mis-type

ALLOWED_FUNCTIONALITIES : set[int] = {2}
BLACKLISTED_MECHANISMS = [
    'imide',
    'vinyl'
]

# PROJECT-SPECIFIC JOB HELPER FUNCTIONS
## JOB LOGGING
def redirect_job_to_logfile(job : Job) -> logging.Logger:
    '''Thin wrapper around redirect_job_to_logfile() which is compatible with Signac jobs'''
    return redirect_to_logfile(
        logfile_path=job.fn(PolymerBuildProject.LOGFILE_NAME),
        logger_name=job.id,
        level=PolymerBuildProject.LOGLEVEL,
        aux_loggers=POLYMERIST_LOGGERS, # make this all loggers?
        # aux_loggers=get_active_loggers(), # get EVERY active logger registered across all Python modules
    )

## JOB HOOKS
def record_operation_start_time(operation_name : str, job : Job) -> None:
    '''Record into the job document when a particular operation began'''
    job.doc.setdefault(PolymerBuildProject.OP_TIME_RECORD_NAME, {})
    job.doc[PolymerBuildProject.OP_TIME_RECORD_NAME].update({f'{operation_name}_start' : time.time()})

def record_operation_duration(operation_name : str, job : Job) -> None:
    '''Record into the job document how long a particular operation took'''
    op_times = job.doc.get(PolymerBuildProject.OP_TIME_RECORD_NAME, {})
    start_time = op_times.pop(f'{operation_name}_start') # look up and withdraw start time
    if start_time is None: # NOTE: dicts in Signac documents are NOT pure Python dicts; their "pop()" method returns None by default and never raises KeyError
        raise ValueError(f'No start time recorded for operation "{operation_name}"; cannot calculate operation duration')
    
    op_duration = time.time() - start_time
    job.doc[PolymerBuildProject.OP_TIME_RECORD_NAME].update({operation_name : op_duration})

## HELPERS FOR OBTAINING OBJECTS BASED ON JOB DATA
def has_nonempty_file(job : Job, filename : str) -> bool:
    '''Check if a job contains a particular file which contains a nonzero amount of information'''
    return job.isfile(filename) and not is_empty(job.fn(filename))

def load_job_rdmol(job : Job, separate_mols : bool=True) -> Chem.Mol:
    '''Helper method for loading an RDKit molecule from the SMILES in a job's statepoint'''
    return sanitized_mol_from_smiles(job.sp.smiles_explicit, separate_mols=separate_mols)

def load_job_rxn(job : Job) -> AnnotatedReaction:
    '''Helper method for loading an RDKit molecule from the SMILES in a job's statepoint'''
    rxn = AnnotatedReaction.from_smarts(job.sp.rxn_smarts.replace('#0', '*'))
    rxn.Initialize()
    n_warn, n_err = rxn.Validate()
    # assert n_err == 0

    return rxn

def load_job_topology(job : Job, sdf_pathname : str, *start_args, **kwargs) -> Optional[Molecule]:
    '''Read and return an OpenFF Topology from well-formed SDF file,
    returning None if encoding or other errors are encountered'''
    if not job.isfile(sdf_pathname):
        return None
    
    with redirect_job_to_logfile(job) as logger:
        sdf_path = Path(job.fn(sdf_pathname))
        if sdf_path.suffix != '.sdf':
            logger.error('Only SDF files are allowed for loading molecules')
            return None

        try:
            return topology.topology_from_sdf(sdf_path, *start_args, **kwargs)
        except UnassignedChemistryInPDBError: # special cases for known common errors
            logger.error('OpenFF will not load molecule with ambiguous stereochemistry')
            return None
        except MoleculeParseError: # special cases for known common errors
            logger.error('Empty or malformed SDF file, could not read structural data')
            return None

load_job_melt_neat_topology = partial(
    load_job_topology,
    sdf_pathname=PolymerBuildProject.MELT_NEAT_SDF,
    allow_undefined_stereo=PolymerBuildProject.RELAXED_STEREO
)
def load_job_oligomer_molecule(job : Job) -> Optional[Molecule]:
    '''Check whether a topology has atomic partial charges assigned to it'''
    oligomer_top = load_job_topology(
        job, PolymerBuildProject.OLIGOMER_SDF,
        allow_undefined_stereo=PolymerBuildProject.RELAXED_STEREO
    )
    if oligomer_top is None:
        return None
    return topology.get_largest_offmol(oligomer_top)

def load_job_interchange(job : Job) -> Optional[Interchange]:
    '''Load a serialized OpenFF Interchange object into memory if one has been written'''
    if not job.isfile(PolymerBuildProject.INTERCHANGE_PATH):
        return None
    
    with open(job.fn(PolymerBuildProject.INTERCHANGE_PATH), 'rb') as file:
        inc = pickle.load(file)
        return inc
    
def load_job_openmm_context(job : Job) -> Context:
    '''Load a complete OpenMM Context from saved System, State, and Integrator'''
    with open(job.fn(PolymerBuildProject.OPENMM_SYSTEM_PATH), 'r') as sys_file, \
        open(job.fn(PolymerBuildProject.OPENMM_STATE_PATH ), 'r') as state_file, \
        open(job.fn(PolymerBuildProject.OPENMM_INTEG_PATH ), 'r') as integ_file:
        system = XmlSerializer.deserialize(sys_file.read())
        state = XmlSerializer.deserialize(state_file.read())
        integrator = XmlSerializer.deserialize(integ_file.read())

    context = Context(system, integrator)
    apply_state_to_context(context, state)

    return context


# OPERATIONS, LABELS AND CONDITIONS
everything = PolymerBuildProject.make_group(name='everything') # "master" group which allows submission of all operations

## 0) CHEMISTRY VALIDATION
validate_chemistry = PolymerBuildProject.make_group(name='validate_chemistry')

def atoms_validated(job : Job) -> bool:
    return 'has_banned_atom_types' in job.doc

def monomer_compositions_validated(job : Job) -> bool:
    return 'has_banned_monomer_compositions' in job.doc

def monomer_sizes_validated(job : Job) -> bool:
    return 'has_oversized_monomers' in job.doc

def chemistry_validated(job : Job) -> bool: # NOTE: this might seem redundant at a glance, but enables delineation between unknown chemical status and KNOWN chemical invalidity
    '''Check if chemistry prechecks have been done (i.e. status is not unknown)'''
    return atoms_validated(job) \
        and monomer_compositions_validated(job) \
        and monomer_sizes_validated(job) \
            
def chemistry_passed(job : Job) -> bool:
    '''Check if all chemical checks have been passed in aggregate'''
    return not ( # chemistry is valid IFF ALL of the below have been evaluated and are NOT found to be invalid
        has_banned_atom_types(job) \
        or has_banned_monomer_compositions(job) \
        or has_oversized_monomers(job) \
    )

### NOTE: found it less cluttered to indicate when chemical conditions are NOT met, rather than dsplay when ALL have passed
@PolymerBuildProject.label 
def has_banned_atom_types(job : Job) -> bool:
    return job.doc.get('has_banned_atom_types', False) # default to False if no substructure evaluation has been performed yet 

@PolymerBuildProject.label
def has_banned_monomer_compositions(job : Job) -> bool:
    return job.doc.get('has_banned_monomer_compositions', False) # default to False if no substructure evaluation has been performed yet 

@PolymerBuildProject.label
def has_oversized_monomers(job : Job) -> bool:
    return job.doc.get('has_oversized_monomers', False) # default to False if no substructure evaluation has been performed yet 

@PolymerBuildProject.label
def chemistry_allowed(job : Job) -> bool:
    '''Check if all chemical checks are KNOWN to be passing'''
    return chemistry_validated(job) and chemistry_passed(job)

@everything
@validate_chemistry
@PolymerBuildProject.post(atoms_validated)
@PolymerBuildProject.operation(directives={'walltime' : 1/60, 'np' : 1})
def validate_atoms(job : Job) -> None:
    '''Check that no banned atoms types are present in any monomer molecules'''
    monomers = load_job_rdmol(job, separate_mols=False) # no need to separate, since checking if ANY submols contain a banned substructure
    with redirect_job_to_logfile(job) as logger:
        logger.info('Searching for banned atom types')
        for query_name, atom_query in PolymerBuildProject.BLACKLISTED_ATOM_QUERIES.items():
            if any(monomers.GetAtomsMatchingQuery(atom_query)):
                job.doc['has_banned_atom_types'] = True
                logger.error(f'Detected invalid atoms of type "{query_name}"') # TODO: make this more descriptive
                break # no need to screen any further queries if even one disallowed result is detected
        else:
            job.doc['has_banned_atom_types'] = False
            logger.info('No banned atom types detected')

@everything
@validate_chemistry
@PolymerBuildProject.pre.not_(has_banned_atom_types) # no point in checking monomer compositions if the more primitive atom type check fails
@PolymerBuildProject.post(monomer_compositions_validated)
@PolymerBuildProject.operation(directives={'walltime' : 1/60, 'np' : 1})
def validate_monomer_compositions(job : Job) -> None:
    '''Check that no banned monomer molecules are present in as any of the monomer input molecules'''
    monomers = load_job_rdmol(job, separate_mols=False) # no need to separate, since checking if ANY submols contain a banned substructure
    with redirect_job_to_logfile(job) as logger:
        logger.info('Searching for banned monomer compositions')
        for smiles, banned_monomer in BLACKLISTED_MONOMER_MOLS.items():
            if monomers.HasSubstructMatch(banned_monomer):
                job.doc['has_banned_monomer_compositions'] = True
                logger.error(f'Detected invalid monomer "{smiles}"') # TODO: make this more descriptive
                break # no need to screen any further queries if even one disallowed result is detected
        else:
            job.doc['has_banned_monomer_compositions'] = False
            logger.info('No banned monomer compositions detected')

@everything
@validate_chemistry
@PolymerBuildProject.pre.not_(has_banned_monomer_compositions) 
@PolymerBuildProject.post(monomer_sizes_validated)
@PolymerBuildProject.operation(directives={'walltime' : 1/60, 'np' : 1})
def validate_monomer_sizes(job : Job) -> None:
    '''Check that none of the monomer input molecules are larger than the prescribed limit'''
    monomers = load_job_rdmol(job, separate_mols=False), # NOTE: here we DO need to separate, since we are looking at the sizes of individual molecules
    with redirect_job_to_logfile(job) as logger:
        logger.info('Searching for oversized monomers')
        for monomer in monomers:
            if (n_atoms := monomer.GetNumAtoms()) > PolymerBuildProject.N_ATOM_CAP_MONOMER:
                logger.error(f'Detected oversized monomer (containing {n_atoms} atoms relative to the prescribed {PolymerBuildProject.N_ATOM_CAP_MONOMER}-atom cutoff)')
                job.doc['has_oversized_monomers'] = True
                break
        else:
            logger.info('No oversized monomers detected')
            job.doc['has_oversized_monomers'] = False
    
        
## 1) REACTANT AND MECHANISM PERCEPTION
perceive_mechanism = PolymerBuildProject.make_group(name='perceive_mechanism')

@PolymerBuildProject.label
def labelled_mechanism_allowed(job : Job) -> bool:
    '''Check that the rxn mechanism type is not explicitly blacklisted'''
    return ('mechanism_labelled' in job.doc) and (job.doc.mechanism_labelled not in BLACKLISTED_MECHANISMS)

def reactant_order_evaluated(job : Job) -> bool:
    return 'reactant_ordering' in job.doc

@PolymerBuildProject.label
def matches_rxn_template(job : Job) -> bool:
    return reactant_order_evaluated(job) and (job.doc.reactant_ordering is not None)

# REACTANT FUNCTIONALITY CHECK
def functionalities_evaluated(job : Job) -> bool:
    return 'functionalities' in job.doc

@PolymerBuildProject.label
def monomers_satisfy_functionality(job : Job) -> bool:
    '''Check that all monomers have allowed degrees of functionalization'''
    return functionalities_evaluated(job) and all(f in ALLOWED_FUNCTIONALITIES for f in job.doc.functionalities)

@PolymerBuildProject.label
def has_chemical_fragments(job : Job) -> bool:
    '''Check if repeat unit fragments from reaction enumeration have been cached'''
    return has_nonempty_file(job, PolymerBuildProject.FRAGMENTS_PATH)

@everything
@perceive_mechanism
@PolymerBuildProject.pre(chemistry_allowed)
@PolymerBuildProject.pre(labelled_mechanism_allowed)
@PolymerBuildProject.post(reactant_order_evaluated)
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def determine_reactant_order(job : Job) -> None:
    '''
    Check that SMILES monomers are compatible with the 
    reaction template for the mechanism they claim to follow
    '''
    # 1) check that monomers fit a reaction template
    reactants = load_job_rdmol(job, separate_mols=True)
    rxn = load_job_rxn(job)

    with redirect_job_to_logfile(job) as logger:
        try:
            reactant_ordering = rxn.valid_reactant_ordering(reactants, as_mols=False)
            job.doc.reactant_ordering = reactant_ordering # set EVEN if found ordering is None, to indicate this check has already been done
        except BadNumberReactants as bnr_error: # temporarily intercept this error to mark the job as having been checked for post-conditions
            job.doc.reactant_ordering = None
            raise bnr_error # re-raise to propagate this error up to the logger context
        
        if reactant_ordering is not None:
            logger.info(f'Identified valid reactant ordering: {reactant_ordering}')
        else:
            logger.error(f'No valid ordering of reactants could be solved for the chosen "{job.doc.mechanism_labelled}" rxn template')

@everything
@perceive_mechanism
@PolymerBuildProject.pre(chemistry_allowed)
@PolymerBuildProject.pre(labelled_mechanism_allowed)
@PolymerBuildProject.pre(matches_rxn_template)
@PolymerBuildProject.post(functionalities_evaluated)
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def determine_reactant_functionalities(job : Job) -> None:
    '''
    Once a reactant ordering has been identified, determine how many of
    each template group are present in each respective reactant monomer
    '''
    rxn = load_job_rxn(job)
    reactant_smiles_all = job.sp.smiles_explicit.split('.')

    with redirect_job_to_logfile(job) as logger:
        functionalities : list[int] = []
        for i in job.doc.reactant_ordering:
            reactant_smiles = reactant_smiles_all[i]
            reactant_mol = sanitized_mol_from_smiles(reactant_smiles)
            
            num_funct_groups = num_substruct_queries_distinct(reactant_mol, rxn.GetReactantTemplate(i))
            functionalities.append(num_funct_groups)
            if num_funct_groups not in ALLOWED_FUNCTIONALITIES:
                logger.error(f'Found {num_funct_groups} active functional groups (vs any from {ALLOWED_FUNCTIONALITIES}) for molecule {reactant_smiles}')

        job.doc.functionalities = functionalities

@everything
@perceive_mechanism
@PolymerBuildProject.pre(matches_rxn_template)
@PolymerBuildProject.pre(monomers_satisfy_functionality)
@PolymerBuildProject.post(has_chemical_fragments)
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def enum_fragments(job : Job) -> None:
    '''Enumerate all possible repeat unit fragment using cheminformatic reaction procedure'''
    reactants = load_job_rdmol(job, separate_mols=True)
    rxn = load_job_rxn(job)
    reactor = PolymerizationReactor(rxn)
    
    monogrp = MonomerGroup()
    with redirect_job_to_logfile(job) as logger:
        monogrp = generate_smarts_fragments(reactants, reactor)
        monogrp.to_file(job.fn(PolymerBuildProject.FRAGMENTS_PATH))
        logger.info('Successfully enumerated and cached repeat unit fragments')


## 2) TOPOLOGY ASSEMBLY AND COORDINATE GENERATION
oligomerize = PolymerBuildProject.make_group(name='oligomerize')

@PolymerBuildProject.label
def coordinates_generated(job : Job) -> bool:
    '''Check whether a nonempty PDB file has been generated from fragments'''
    return has_nonempty_file(job, PolymerBuildProject.OLIGOMER_PDB)

@PolymerBuildProject.label
def matches_m2p_smiles(job : Job) -> bool:
    '''Check whether the resulting polymer agrees with the SMILES output of M2P (if data is provided)'''
    ... # this check is not always applicatble, and the way DOP is currently defined is a little bit dicey

@PolymerBuildProject.label
def no_ring_piercing(job : Job) -> bool:
    ... # TODO: implement post-minimization bond length check
    
@everything
@oligomerize
# @PolymerBuildProject.pre.copy_from(determine_reactant_order)
# @PolymerBuildProject.pre.copy_from(enum_fragments)
@PolymerBuildProject.pre(has_chemical_fragments)
@PolymerBuildProject.post(coordinates_generated)
@PolymerBuildProject.operation(directives={'walltime' : 15/60, 'np' : 1})
def build_oligomer_pdb(job : Job) -> None:
    '''Generate coordinates and build oligomer PDB file using mBuild'''
    monogrp = MonomerGroup.from_file(job.fn(PolymerBuildProject.FRAGMENTS_PATH))
    with redirect_job_to_logfile(job) as logger:
        # check for identical parallel oligomer jobs
        parallel_struct_jobs = job.project.find_jobs({
            f'sp.{attr}' : getattr(job.sp, attr)
                for attr in ('DOP', 'smiles_explicit')
        })
        for parallel_job in parallel_struct_jobs:
            if (parallel_job.id != job.id) and (parallel_job.isfile(PolymerBuildProject.OLIGOMER_PDB)):
                logger.info(f'Job {parallel_job.id} already has my oligomer PDB!')
                break

        # generate coordinates with mBuild hook
        polymer = build_linear_polymer(
            monomers=monogrp,
            n_monomers=(job.sp.DOP*len(job.sp.copolymer_sequence)), # interpret DOP here as number of monomer sequence repeats (including end groups)
            sequence=job.sp.copolymer_sequence, # DEV: for now, fixed as "BA" for all mechanisms; TODO: find way to set this as a function of mechanism in setup
            energy_minimize=True, # TODO: add master config option for energy minimization at project level
        )
        mbmol_to_openmm_pdb(job.fn(PolymerBuildProject.OLIGOMER_PDB), polymer)
        logger.info('Successfully generated PDB structure file')

### 2A) OPENFF PARAMETER ASSIGNMENT
@PolymerBuildProject.label
def chemical_info_assigned(job : Job) -> bool:
    '''Check whether a topology has atomic partial charges assigned to it'''
    return has_nonempty_file(job, PolymerBuildProject.OLIGOMER_SDF)# and is_valid_sdfile(job.fn(PolymerBuildProject.OLIGOMER_SDF)) #(load_job_oligomer_molecule(job) is not None)

@PolymerBuildProject.label
def partial_charges_assigned(job : Job) -> bool:
    '''Check whether a topology has atomic partial charges assigned to it'''
    return 'partial_charges' in job.data

@everything
@oligomerize
# @PolymerBuildProject.pre.copy_from(build_oligomer_pdb)
@PolymerBuildProject.pre(coordinates_generated)
@PolymerBuildProject.post(chemical_info_assigned)
@PolymerBuildProject.operation(directives={'walltime' : 20/60, 'np' : 1})
def assign_chem_info(job : Job) -> None:
    '''Assign chemical information to bare PDB graph and export completely-specified system to SDF file'''
    monogrp = MonomerGroup.from_file(job.fn(PolymerBuildProject.FRAGMENTS_PATH))
    with redirect_job_to_logfile(job) as logger:
        offtop = Topology.from_pdb(job.fn(PolymerBuildProject.OLIGOMER_PDB), _custom_substructures=monogrp.monomers)
        if not partition(offtop):
            logger.error(f'Failed to produce residue partition with fragments for job {job.id}')
            return None # exit before writing SDF; will cause post-condition to not be meet
        topology.topology_to_sdf(job.fn(PolymerBuildProject.OLIGOMER_SDF), offtop)
        logger.info('Successfully generated chemically-explicit SDF structure file')

@everything
@oligomerize
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.post.true('r_eff') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('n_atoms_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('elem_counts_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('molar_mass_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def summarize_oligomer(job : Job) -> None:
    '''Compute simple summarizing info about an oligomer which streamline lattice packing,
    namely the number of atoms, the distribution of elements and the effect (max) radius'''
    offmol = load_job_oligomer_molecule(job)
    job.doc.r_eff = rdmol_effective_radius(offmol.to_rdkit()) # TODO: this shouldn't require an RDKit conversion, just access tot he conformer
    job.doc.n_atoms_oligomer = offmol.n_atoms
    job.doc.elem_counts_oligomer = elem_counts(offmol)
    job.doc.molar_mass_oligomer = sum(atom.mass for atom in offmol.atoms).magnitude
    
@everything
@oligomerize
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.post(partial_charges_assigned) # TODO: fill this in with something more substantive!!
@PolymerBuildProject.operation(directives={'walltime' : 3/60, 'np' : 1})
def assign_partial_charges(job : Job) -> None:
    '''Generate coordinates and build oligomer PDB file using mBuild'''
    offmol = load_job_oligomer_molecule(job)
    charger_type = MolCharger.subclass_registry.get(job.sp.pcharge_method, None)
    if charger_type is None:
        return
        
    with redirect_job_to_logfile(job) as logger:
        charger = charger_type()
        cmol = charger.charge_molecule(offmol)
        topology.topology_to_sdf(job.fn(PolymerBuildProject.OLIGOMER_SDF), cmol.to_topology())

        # cache partial charge data to job records
        pcharge_unit : OFFUnit = cmol.partial_charges.units
        job.data.partial_charges = cmol.partial_charges.m_as(pcharge_unit)
        job.doc.pcharge_units = f'{pcharge_unit:simple}' # convert to string with explicit formatting to allow recovery of Unit type from text


## 3) LATTICE SIZING AND PACKING
pack_lattice = PolymerBuildProject.make_group(name='pack_lattice') 

def lattice_sites_determined(job : Job) -> bool:
    '''Check whether lattice sites have been assigned for '''
    return 'lattice_sites' in job.data

@PolymerBuildProject.label
def neat_melt_packed(job : Job) -> bool:
    '''Check is packing of the neat melt was successful'''
    return has_nonempty_file(job, PolymerBuildProject.MELT_NEAT_SDF)# and is_valid_sdfile(job.fn(PolymerBuildProject.MELT_NEAT_SDF)) #and (load_job_melt_neat_topology(job) is not None)

def pbcs_determined(job : Job) -> bool:
    '''Check whether periodic bounding box has been calculated'''
    return ('box_vectors_nm' in job.data) and ('box_vector_dims' in job.doc)

@everything
@pack_lattice
@PolymerBuildProject.pre(chemical_info_assigned) # don't need charges, only valid cornformer to pick sites
@PolymerBuildProject.post.true('n_oligomers')
@PolymerBuildProject.post.true('lattice_shape')
@PolymerBuildProject.post(lattice_sites_determined)
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def determine_lattice_sites(job : Job) -> None:
    '''Choose smallest accomodating cubic lattice, randomly subsample sites, and scale appropriately to oligomer size'''
    int_lattice = generate_uniform_subpopulated_lattice(
        max_num_atoms=job.sp.n_atoms_max,
        num_atoms_in_mol=job.doc.n_atoms_oligomer
    )
    job.doc.n_oligomers = int_lattice.n_points
    job.doc.lattice_shape = int_lattice.counts_along_dims_as_str()
    
    transform = 2.0 * job.doc.r_eff * np.eye(3, dtype=float) # uniform scaling of lattice which guarantees points are one effective diameter apart
    job.data.lattice_sites = int_lattice.linear_transformation(transform, as_coords=False) # save lattice sites to numpy array on disc

@everything
@pack_lattice
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.pre(partial_charges_assigned)
@PolymerBuildProject.pre(lattice_sites_determined)
@PolymerBuildProject.post(neat_melt_packed)
@PolymerBuildProject.operation(directives={'walltime' : 15/60, 'np' : 1})
def pack_oligomers_onto_lattice(job : Job) -> None:
    '''Clone, randomly rotate, and move oligomer onto predetermined lattice sites'''
    offmol = load_job_oligomer_molecule(job)
    with job.data:
        lattice_sites = job.data.lattice_sites[:]

    with redirect_job_to_logfile(job) as logger:
        melt_offtop = topology.topology_from_molecule_onto_lattice(
            offmol,
            lattice_points=lattice_sites,
            rotate_randomly=True,
            unique_mol_ids=True
        )
        topology.topology_to_sdf(job.fn(PolymerBuildProject.MELT_NEAT_SDF), melt_offtop)

@everything
@pack_lattice
@PolymerBuildProject.pre(neat_melt_packed)
@PolymerBuildProject.post(pbcs_determined)
@PolymerBuildProject.operation(directives={'walltime' : 2/60, 'np' : 1})
def determine_periodic_box(job : Job) -> None:
    '''Size and cache periodic box vectors for the resulting melt topology'''
    melt_offtop    : Topology    = load_job_melt_neat_topology(job)
    box_padding    : OFFQuantity = job.sp.box_padding_nm * offunit.nanometer
    nonbond_cutoff : OFFQuantity = job.sp.nonbonded_cutoff_nm * offunit.nanometer
    min_bbox       : OFFQuantity = 2 * nonbond_cutoff * np.eye(3) # box should be at least twice the nonbonded cutoff to avoid self-interaction

    with redirect_job_to_logfile(job) as logger:
        melt_box_vectors = boxvectors.get_topology_bbox(melt_offtop) # determine tight box size
        melt_box_vectors = boxvectors.pad_box_vectors_uniform(melt_box_vectors, box_padding) # pad out by target amount
        melt_box_vectors = np.maximum(min_bbox, melt_box_vectors) # ensure box is no smaller than the minimum allowed size
        if isinstance(melt_box_vectors, OMMQuantity):
            melt_box_vectors = openmm_to_openff(melt_box_vectors) # ensforce Pint units to to nice m_as conversion

        # cache box vectors and summary
        box_vector_dims : dict[str, str] = {}
        box_vector_str_parts : list[str] = []

        box_vector_unit : OFFUnit = melt_box_vectors.units
        for axis, box_dim_quantity in zip('xyz', np.linalg.norm(melt_box_vectors, axis=1)):
            box_dim_value = round(box_dim_quantity.m_as(box_vector_unit), PolymerBuildProject.QUANTITY_PRECISION)
            box_vector_dims[f'Box dim {axis}'] = f'{(box_dim_value*box_vector_unit)!s}'
            box_vector_str_parts.append(str(box_dim_value))

        box_dim_summary_str = ' x '.join(box_vector_str_parts) + f' {box_vector_unit!s}'
        logger.info(f'Setting {box_dim_summary_str} periodic box')

        job.doc.box_vector_dims = box_vector_dims
        job.data.box_vectors_nm = melt_box_vectors.m_as(offunit.nanometer) # store just the array of vectors (no units) in nm


## 4) OPENFF INTERCHANGE EXPORT
to_interchange = PolymerBuildProject.make_group(name='to_interchange') 

@PolymerBuildProject.label
def forcefield_is_valid(job : Job) -> bool:
    '''Check that the force field specified is a valid and loadable OpenFF forcefield file installed in the current environment'''
    try: # TODO: worth checking explicitly that the file exists/sanitizing missing .offxml etc.?
        ForceField(job.sp.forcefield) # NOTE: need to handle exception when the offxml provided doesn't exist
        return True
    except OSError as error: # TODO: make error handling more specific and informative, left suggestive of common OSError for now
        return False
    
@PolymerBuildProject.label
def has_interchange(job : Job) -> bool:
    '''Whether the current job has produced and serialized an OpenFF Interchange object'''
    return has_nonempty_file(job, PolymerBuildProject.INTERCHANGE_PATH)

@PolymerBuildProject.label # NOTE: the way this label is set up is potentially confusion, but is motivatated by the fact that we only want to flag when stereo does NOT agreed (i.e. be silent when it does/has not been determined)
def interchange_stereo_inconsistent(job : Job) -> bool:
    '''Indicate that interchange creation could not be performed due to library charge mismatch from stereo incompatibility (a known bug)'''
    return not job.doc.get('interchange_stereo_consistent', True) # assumes consistency when it has not yet been explicitly determined (only returns True overall is explicitly found to be inconsistent)

@everything
@to_interchange
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.pre(partial_charges_assigned)
@PolymerBuildProject.pre(neat_melt_packed)
@PolymerBuildProject.pre(pbcs_determined)
@PolymerBuildProject.pre(forcefield_is_valid)
@PolymerBuildProject.pre.not_(interchange_stereo_inconsistent)
@PolymerBuildProject.post(has_interchange)
@PolymerBuildProject.operation(directives={'walltime' : 20/60, 'np' : 1})
def neat_melt_to_interchange(job : Job) -> None:
    '''Create Interchange from final melt (w/ appropriate FF parameters and cutoffs) and pickle for reuse'''
    cmol = load_job_oligomer_molecule(job) # need charged molecule for reference to avoid expensive AM1-BCC default
    melt_offtop = load_job_melt_neat_topology(job)

    with redirect_job_to_logfile(job) as logger:
        nonbond_cutoff = job.sp.nonbonded_cutoff_nm * offunit.nanometer
        logger.info(f'Using {nonbond_cutoff!s} nonbonded interaction cutoff')
        with job.data:
            box_vectors = job.data.box_vectors_nm[:] * offunit.nanometer

        if job.sp.use_switching_function:
            switch_width = job.sp.switch_width_nm*offunit.nanometer
            logger.info(f'Enabling nonbonded switching function with cutoff distance of {switch_width!s}')
        else:
            switch_width = 0.0*offunit.nanometer
            logger.warning('Disabling switching function for nonbonded forces')

        logger.info(f'Obtaining force field parameters from OpenFF "{job.sp.forcefield}"')
        forcefield = ForceField(job.sp.forcefield) # NOTE: need to handle exception when the offxml provided doesn't exist
        if 'ToolkitAM1BCC' in forcefield.registered_parameter_handlers:
            forcefield.deregister_parameter_handler('ToolkitAM1BCC') # forcibly remove AM1BCC handler so a fail siomorphism doesn't result in prohibitively-long AM1BCC calculation

        try:
            logger.info('Initializing OpenFF Interchange from melt topology and force field')
            interchange = forcefield.create_interchange(melt_offtop, charge_from_molecules=[cmol])
        except RuntimeError:
            logger.error('Could not create OpenFF Interchange instance due to stereochemistry incompatibility with single oligomer')
            job.doc.interchange_stereo_consistent = False
            return
        else:
            logger.info('Successfully created OpenFF Interchange instance')
            job.doc.interchange_stereo_consistent = True
        
        # configure PBCs and nonbonded parameters
        interchange.box = box_vectors
        interchange['vdW'].switch_width = switch_width
        interchange['vdW'].cutoff = nonbond_cutoff
        interchange['Electrostatics'].cutoff = nonbond_cutoff

        # pickle for reuse
        with open(job.fn(PolymerBuildProject.INTERCHANGE_PATH), 'wb') as pklfile: # NOTE: pickled files must be read/written in binary mode
            logger.info('Pickling Interchange for reuse in MD export')
            pickle.dump(interchange, pklfile)


## 5) EXPORT TO MD ENGINES
md_export = PolymerBuildProject.make_group(name='md_export') 
openmm_export = PolymerBuildProject.make_group(name='openmm_export') 
lammps_export = PolymerBuildProject.make_group(name='lammps_export') 

### LAMMPS
@PolymerBuildProject.label
def exported_to_lammps(job : Job) -> bool:
    '''Check if LAMMPS files have been generated'''
    return all(
        has_nonempty_file(job, lmp_file)
            for lmp_file in PolymerBuildProject.LAMMPS_PATHS
    )

@PolymerBuildProject.label
def energies_evaluated_lammps(job : Job) -> bool:
    '''Check if LAMMPS energy evaluation has been performed'''
    return has_nonempty_file(job, PolymerBuildProject.LAMMPS_ENERGIES)

@everything
@md_export
@lammps_export
@PolymerBuildProject.pre(has_interchange)
@PolymerBuildProject.post(exported_to_lammps)
@PolymerBuildProject.operation(directives={'walltime' : 10/60, 'np' : 1})
def export_lammps_files(job : Job) -> None:
    '''Write LAMMPS data and input files for the melt topology'''
    interchange = load_job_interchange(job)
    Path(job.fn(PolymerBuildProject.LAMMPS_DIR)).mkdir(exist_ok=True) # ensure the child directory exists

    with redirect_job_to_logfile(job) as logger:
        logger.info('Writing LAMMPS data file')
        interchange.to_lammps_datafile(job.fn(PolymerBuildProject.LAMMPS_DATA_PATH))
        
        logger.info('Writing LAMMPS input file')
        interchange.to_lammps_input(
            file_path=job.fn(PolymerBuildProject.LAMMPS_INPUT_PATH),
            data_file=f'"./{PolymerBuildProject.LAMMPS_DATA_PATH}"', # write reference to data file relative to workspace directory
        )
        logger.info('LAMMPS files successfully written')

@everything
@md_export
@lammps_export
@PolymerBuildProject.pre(exported_to_lammps)
@PolymerBuildProject.post(energies_evaluated_lammps)
@PolymerBuildProject.operation(directives={'walltime' : 5/60, 'np' : 1})
def evaluate_energies_LAMMPS(job : Job) -> None:
    '''Evaluate starting structure energies of LAMMPS MD files'''
    LMP_ARGS = ["-screen", "none", "-log", "none"] # blocks stdout and log.lammps writes to avoid clutter
    with redirect_job_to_logfile(job) as logger:
        with job: # need this context to make relative path work
            energies_lmp_raw = get_lammps_energies(
                job.fn(PolymerBuildProject.LAMMPS_INPUT_PATH),
                preferred_unit=PolymerBuildProject.ENERGY_UNIT,
                cmdargs=LMP_ARGS,
            )
        logger.info('Completed energy evaluation from LAMMPS input files')

        energies_lmp = { # reformat and combine into standard form
            'Potential' : energies_lmp_raw['Potential'],
            'Bond'      : energies_lmp_raw['Bond'],
            'Angle'     : energies_lmp_raw['Angle'],
            'Dihedral'  : energies_lmp_raw['Proper Torsion'] + energies_lmp_raw['Improper Torsion'],
            'vdW'       : energies_lmp_raw['vdW'] + energies_lmp_raw['Dispersion'],
            'Coulomb'   : energies_lmp_raw['Coulomb Short'] + energies_lmp_raw['Coulomb Long'],
        }
        energies_lmp_stringy = {e_name : f'{e_val!s}' for e_name, e_val in energies_lmp.items()} # stringify for JSON serialization
        logger.info('Reformatted LAMMPS energy output')
        
        with open(job.fn(PolymerBuildProject.LAMMPS_ENERGIES), 'w') as energies_lmp_file:
            json.dump(energies_lmp_stringy, energies_lmp_file, indent=4)
        logger.info('LAMMPS energy output saved to file')

### OpenMM
@PolymerBuildProject.label
def exported_to_openmm(job : Job) -> bool:
    '''Check if OpenMM files have been generated'''
    return all(
        has_nonempty_file(job, lmp_file)
            for lmp_file in PolymerBuildProject.OPENMM_PATHS
    )

@PolymerBuildProject.label
def energies_evaluated_openmm(job : Job) -> bool:
    '''Check if OpenMM energy evaluation has been performed'''
    return has_nonempty_file(job, PolymerBuildProject.OPENMM_ENERGIES)

# @everything
@md_export
@openmm_export
@PolymerBuildProject.pre(has_interchange)
@PolymerBuildProject.post(exported_to_openmm)
@PolymerBuildProject.operation(directives={'walltime' : 15/60, 'np' : 1})
def export_openmm_files(job : Job) -> None:
    '''Write OpenMM data and input files for the melt topology'''
    interchange = load_job_interchange(job)
    Path(job.fn(PolymerBuildProject.OPENMM_DIR)).mkdir(exist_ok=True) # ensure the child directory exists

    with redirect_job_to_logfile(job) as logger:
        integrator = LangevinMiddleIntegrator(300*kelvin, 1*picosecond**-1, 1*femtosecond) # hard-coded Intergrator now, will add support for more targetted integrator later
        logger.info('Writing OpenMM topology, system, state, and integrator files')
        interchange_to_openmm(
            interchange=interchange,
            integrator=integrator,
            omm_top_path=job.fn(PolymerBuildProject.OPENMM_TOPO_PATH),
            omm_sys_path=job.fn(PolymerBuildProject.OPENMM_SYSTEM_PATH),
            omm_state_path=job.fn(PolymerBuildProject.OPENMM_STATE_PATH),
            omm_integ_path=job.fn(PolymerBuildProject.OPENMM_INTEG_PATH),
        )
        logger.info('OpenMM files successfully written')

# @everything
@md_export
@openmm_export
@PolymerBuildProject.pre(exported_to_lammps)
@PolymerBuildProject.post(energies_evaluated_openmm)
@PolymerBuildProject.operation(directives={'walltime' : 5/60, 'np' : 1})
def evaluate_energies_openmm(job : Job) -> None:
    '''Evaluate starting structure energies of OpenMM MD files'''
    context = load_job_openmm_context(job)
    with redirect_job_to_logfile(job) as logger:
        omm_pot, omm_kin = eval_openmm_energies_separated(context, preferred_unit=PolymerBuildProject.ENERGY_UNIT)
        logger.info('Completed energy evaluation from OpenMM Context')
        
        KE_total = omm_kin['Total kinetic energy']
        if KE_total != 0.0*PolymerBuildProject.ENERGY_UNIT:
            logger.error(f'Unexpected nonzero KE evaluated ({KE_total!s})')
            return

        omm_pot_raw = {
            e_name.removesuffix(' potential energy').removesuffix(' force') : e_val
                for e_name, e_val in omm_pot.items()
        }
        energies_omm = {
            'Potential' : omm_pot_raw['Total'],
            'Bond'      : omm_pot_raw['HarmonicBondForce'],
            'Angle'     : omm_pot_raw['HarmonicAngleForce'],
            'Dihedral'  : omm_pot_raw['PeriodicTorsionForce'],
            'vdW'       : omm_pot_raw['vdW'] + omm_pot_raw['vdW 1-4'],
            'Coulomb'   : omm_pot_raw['Electrostatics'] + omm_pot_raw['Electrostatics 1-4'],
        }
        energies_omm_stringy = {e_name : f'{e_val!s}' for e_name, e_val in energies_omm.items()} # stringify for JSON serialization
        logger.info('Reformatted OpenMM energy output')
        
        with open(job.fn(PolymerBuildProject.OPENMM_ENERGIES), 'w') as energies_omm_file:
            json.dump(energies_omm_stringy, energies_omm_file, indent=4)
        logger.info('OpenMM energy output saved to file')


# enabling CLI interaction
def main() -> None:
    parser = ArgumentParser()
    parser.add_argument(
        '-path',
        '--project-path',
        type=Path,
        default=Path.cwd(),
        required=True,
        help='Path to the directory in which the (presumed initialized) Signac project statepoints reside',
    ),
    parser.add_argument( 
        '-rxns',
        '--rxn-mapping-path',
        type=Path,
        default=RXNS_DIR/'rxns_polyID.json',
        help='The path to a JSON file containing a reaction name mapping\n' \
            'Should contain dict whose keys are reaction names and whose values are reaction SMARTS strings'
    )
    parser.add_argument(
        '-qp',
        '--quantity-precision',
        type=int,
        default=4,
        help='The number of decimal places to which to display and log physical and numeric quantities',
    )
    parser.add_argument(
        '-arom',
        '--aromaticity-model',
        choices=AromaticityModel.names.keys(),
        default='AROMATICITY_MDL',
        help='The aromaticity model to use when perceiving aromatic bonds',
    )
    parser.add_argument(
        '-sanops',
        '--sanitization-operations',
        choices=SanitizeFlags.names.keys(),
        default='SANITIZE_ALL',
        help='The chemical cleanup operations to be performed any time a molecule is loaded from SMILES',
    )
    parser.add_argument(
        '--strict-stereo',
        action='store_true',
        help='Optional, whether to enforce strict and unambiguous stereochemistry when loading molecules into OpenFF toolkit objects',
    )
    parser.add_argument(
        '-natmcap',
        '--n-atom-cap',
        type=int,
        default=150,
        help='A cap on the number of atoms any individual monomer molecule contains; any chemistries with monomers larger than this cap will NOT be built!'
    )
    # TODO : implement log level setting

    # separate this script's args from those required by signac
    start_args, signac_args = parser.parse_known_args()
    assert start_args.project_path.exists() and start_args.project_path.is_dir()
    
    # configure global vars in Project definition and initialize project instance
    PolymerBuildProject.QUANTITY_PRECISION = start_args.quantity_precision
    PolymerBuildProject.RELAXED_STEREO = not start_args.strict_stereo
    PolymerBuildProject.REGISTERED_RXN_SMARTS = read_rxn_mapping_data(start_args.rxn_mapping_path)
    
    PolymerBuildProject.N_ATOM_CAP_MONOMER = start_args.n_atom_cap
    ## NOTE: these will raise a KeyError if an invalid model name is provided
    PolymerBuildProject.SANITIZE_OPS = SanitizeFlags.names[start_args.sanitization_operations]
    PolymerBuildProject.AROMATICITY_MODEL = AromaticityModel.names[start_args.aromaticity_model]
    
    # PolymerBuildProject.LOGLEVEL = ...
    logging.basicConfig(level=PolymerBuildProject.LOGLEVEL, force=True)

    new_project = PolymerBuildProject(
        path=start_args.project_path,
        entrypoint={
            'path' : __file__,
        }
    )
    # register project-level hooks
    new_project.project_hooks.on_start = [record_operation_start_time]
    new_project.project_hooks.on_exit  = [record_operation_duration  ]
    # new_project.project_hooks.on_exception = ... # TODO: add distinct timing mode for cases where operations fail before exit

    # mock remaining Signac args for parser and run Project's shell interface
    sys.argv[1:] = signac_args # NOTE: this is an ugly hack to allow this script to take CLI args while not disturbing Signacs tastes for arguments
    new_project.main()

if __name__ == '__main__':
    main()
