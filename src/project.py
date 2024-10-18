'''
The structure of a polymer building project workflow,
including a shared namespace, labels, conditions, and operations
'''

import sys
import logging
from argparse import ArgumentParser, Namespace, REMAINDER

from typing import ClassVar, Optional

import pickle
from functools import partial
from pathlib import Path
from string import ascii_uppercase

import numpy as np

# OpenFF toolkits
from openff.toolkit import Molecule, Topology, ForceField
from openff.toolkit.utils.exceptions import (
    MoleculeParseError,
    UnassignedChemistryInPDBError,
    IncorrectNumConformersWarning,
    InconsistentStereochemistryError,
)
logging.getLogger('openff.toolkit.typing.engines.smirnoff.parameters').setLevel(logging.CRITICAL) # silence annoying Electrostatics up-conversion INFO logs

from openff.units import (
    unit as offunit,
    Unit as OFFUnit,
    Quantity as OFFQuantity,
)

from openff.interchange import Interchange

# Cheminformatics
from rdkit import Chem

# OpenMM imports
from openmm.unit import Quantity as OMMQuantity
from openmm import VerletIntegrator
from openmm.unit import femtosecond

# Signac imports
from signac.job import Job
from flow import FlowProject

# Custom (polymerist) imports
import polymerist as ps
from polymerist.genutils.importutils import submodule_loggers
POLYMERIST_LOGGERS = [logger for logger in submodule_loggers(ps).values() if logger is not None] # TODO: move this into polymerist?

from polymerist.unitutils.interop import openmm_to_openff
from polymerist.smileslib import substructures

from polymerist.polymers.monomers import MonomerGroup, specification
from polymerist.polymers.building import build_linear_polymer, mbmol_to_openmm_pdb

from polymerist.mdtools.openfftools import topology, boxvectors
from polymerist.mdtools.openfftools.partition import partition
from polymerist.mdtools.openfftools.partialcharge.molchargers import MolCharger

from polymerist.rdutils.rdcoords.tiling import rdmol_effective_radius
from polymerist.rdutils.reactions.reactions import AnnotatedReaction, BadNumberReactants
from polymerist.rdutils.reactions.reactors import PolymerizationReactor

# Utils imports - made these non-relative to avoid screwing up external vs internal call
from .utils.logs import redirect_to_logfile
from .utils.filelib import is_empty
from .utils.offlib import elem_counts
from .utils.cheminf import is_valid_sdfile
from .utils.packing import generate_uniform_subpopulated_lattice
from .utils.mdexport import interchange_to_lammps, interchange_to_openmm


# ATOMS, MONOMERS, AND REACTION MECHANISMS WHICH ARE, FOR ONE REASON OR ANOTHER, NOT ALLOWED
ALLOWED_FUNCTIONALITIES : set[int] = {2}
BLACKLISTED_ATOM_QUERIES = {
    'sulfur'  : Chem.MolFromSmarts('[S]'),
    # 'phosphorus' : Chem.MolFromSmarts('[P]'),
    'silicon' : Chem.MolFromSmarts('[Si]'),
    'metal'   : substructures.SPECIAL_QUERY_MOLS['metal'],
    # 'halogen' : substructures.SPECIAL_QUERY_MOLS['halogen'],
}
_blacklisted_monomer_smiles = [ # monomers which are, for one reason or another, disallowed
    'CC(C)(C)c1cc(c(Oc2ccc(cc2)N(c3ccc(N)cc3)c4ccc(N)cc4)c(c1)C(C)(C)C)C(C)(C)C',  # the extraordinary number of symmetries of this amine ("4-N-(4-aminophenyl)-4-N-[4-(2,4,6-tritert-butylphenoxy)phenyl]benzene-1,4-diamine")... 
    'CC(C)(C)c1cc(Oc2ccc(-c3ccc(N)cc3)cc2C(F)(F)F)c(C(C)(C)C)cc1Oc1ccc(-c2ccc(N)cc2)cc1C(F)(F)F', # ...mean it takes impractically long to isomorphism match during the Topology partition step
    'CCCCCCCCCCCCCCCCC(CO)C(CO)CCCCCCCCCCCCCCCC', # this one is not necessarily highly-automorphic, but DOES hang up the partition algorithm
] # TODO: might try setting limit of <1000 automorphisms for automatic check (since this is the default limit for substructure matches)
BLACKLISTED_MONOMER_QUERIES = {}
for smiles in _blacklisted_monomer_smiles:
    exp_spi = specification.expanded_SMILES(smiles, assign_map_nums=False)
    banned_mol = Chem.MolFromSmiles(exp_spi, sanitize=False)
    BLACKLISTED_MONOMER_QUERIES[smiles] = banned_mol
BLACKLISTED_MECHANISMS = [
    'imide',
    'vinyl'
]


# DEFINING THE SIGNAC PROJECT CLASS PROPER 
class PolymerBuildProject(FlowProject):
    '''Project for automated high-throughput generation of polymer structure and MD inputs from chemical data'''
    # GLOBAL CONFIG - TODO: make these configuratble via argparse to the containing script
    QUANTITY_PRECISION  : ClassVar[int] = 4 # number of decimal place to report Quantities to when printing/logging
    RELAXED_STEREO      : ClassVar[bool] = True
    LOGLEVEL            : int = logging.INFO

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

    ## OpenMM
    OPENMM_DIR          : ClassVar[str] = 'OpenMM'
    OPENMM_STATE_PATH   : ClassVar[str] = f'{OPENMM_DIR}/state.xml'
    OPENMM_SYSTEM_PATH  : ClassVar[str] = f'{OPENMM_DIR}/system.xml'
    OPENMM_TOPO_PATH    : ClassVar[str] = f'{OPENMM_DIR}/topology.xml'
    OPENMM_INTEG_PATH   : ClassVar[str] = f'{OPENMM_DIR}/integrator.xml'
    OPENMM_PATHS : ClassVar[list[str]] = (
        OPENMM_STATE_PATH,
        OPENMM_SYSTEM_PATH,
        OPENMM_TOPO_PATH,
        OPENMM_INTEG_PATH,
    )


# PROJECT SPECIFIC JOB HELPER FUNCTIONS
def redirect_job_to_logfile(job : Job) -> logging.Logger:
    '''Thin wrapper around redirect_job_to_logfile() which is compatible with Signac jobs'''
    return redirect_to_logfile(
        logfile_path=job.fn(PolymerBuildProject.LOGFILE_NAME),
        logger_name=job.id,
        level=PolymerBuildProject.LOGLEVEL,
        aux_loggers=POLYMERIST_LOGGERS # make this all loggers?
    )

def has_nonempty_file(job : Job, filename : str) -> bool:
    '''Check if a job contains a particular file which contains a nonzero amount of information'''
    return job.isfile(filename) and not is_empty(job.fn(filename))

def load_job_rdmol(job : Job) -> Chem.Mol:
    '''Helper method for loading an RDKit molecule from the SMILES in a job's statepoint'''
    reactant_mol = Chem.MolFromSmiles(job.sp.smiles_explicit, sanitize=False) # CRITICAL that sanitize=False to avoid stripping
    Chem.SanitizeMol(reactant_mol, sanitizeOps=specification.SANITIZE_AS_KEKULE) # single, unified mol containing individual reactant as disconnected components

    return reactant_mol

def load_job_rxn(job : Job) -> AnnotatedReaction:
    '''Helper method for loading an RDKit molecule from the SMILES in a job's statepoint'''
    rxn = AnnotatedReaction.from_smarts(
        job.sp.rxn_smarts.replace('#0', '*') # NOTE: this is a hack which should be sanitized in the previous rxn assembly step
    )
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


# LABELS AND CONDITIONS
## CHEMISTRY PRECHECKS
def atoms_allowed(job : Job) -> bool:
    '''Check no banned atom types are present'''
    # 2) check that none of the monomers are blacklisted
    reactant_mol = load_job_rdmol(job)
    return not any(
        substructures.matching_labels_from_substruct_dict(
            reactant_mol,
            BLACKLISTED_ATOM_QUERIES,
        )
    ) # if any illegal atoms are detected in the current monomer, return and exit

def monomers_allowed(job : Job) -> bool:
    '''Check no banned atom monomer fragments are present'''
    reactant_mol = load_job_rdmol(job)
    return not any(
        substructures.matching_labels_from_substruct_dict(
            reactant_mol,
            BLACKLISTED_MONOMER_QUERIES
        )
    ) # Exclude any monomers which are structurally disallowed

def mechanism_provided(job : Job) -> bool:
    '''Check that a reference "mechanism" field was generated during job initialization'''
    return 'mechanism' in job.doc

def mechanism_allowed(job : Job) -> bool:
    '''Check that the rxn mechanism type is not explicitly blacklisted'''
    return job.doc.mechanism not in BLACKLISTED_MECHANISMS

@PolymerBuildProject.label
def chemistry_valid(job : Job) -> bool:
    '''Aggregate together all atom, monomer, and mechanism prechecks'''
    return atoms_allowed(job) \
        and monomers_allowed(job) \
        and mechanism_provided(job) \
        and mechanism_allowed(job) \

## REACTANT ORDER PERCEPTION
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
    
## STRUCTURE BUILDING VALIDATION
@PolymerBuildProject.label
def has_chemical_fragments(job : Job) -> bool:
    '''Check if repeat unit fragments from reaction enumeration have been cached'''
    return has_nonempty_file(job, PolymerBuildProject.FRAGMENTS_PATH)

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

## OPENFF PARAMETER ASSIGNMENT
@PolymerBuildProject.label
def chemical_info_assigned(job : Job) -> bool:
    '''Check whether a topology has atomic partial charges assigned to it'''
    return has_nonempty_file(job, PolymerBuildProject.OLIGOMER_SDF) and is_valid_sdfile(job.fn(PolymerBuildProject.OLIGOMER_SDF)) #(load_job_oligomer_molecule(job) is not None)

@PolymerBuildProject.label
def partial_charges_assigned(job : Job) -> bool:
    '''Check whether a topology has atomic partial charges assigned to it'''
    if not chemical_info_assigned(job):
        return False
    
    offmol = load_job_oligomer_molecule(job)
    if offmol is None:
        return False # TOSELF: the above check may make this redundant - further testing needed
    
    return offmol.partial_charges is not None

## LATTICE SIZING AND PACKING
def lattice_sites_determined(job : Job) -> bool:
    '''Check whether lattice sites have been assigned for '''
    return 'lattice_sites' in job.data

@PolymerBuildProject.label
def neat_melt_packed(job : Job) -> bool:
    '''Check is packing of the neat melt was successful'''
    return has_nonempty_file(job, PolymerBuildProject.MELT_NEAT_SDF) and is_valid_sdfile(job.fn(PolymerBuildProject.MELT_NEAT_SDF)) #and (load_job_melt_neat_topology(job) is not None)

def pbcs_determined(job : Job) -> bool:
    '''Check whether periodic bounding box has been calculated'''
    return ('box_vectors_nm' in job.data) and ('box_vector_dims' in job.doc)

## OPENFF INTERCHANGE EXPORT
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


# OPERATIONS
# 1) TEST FOR RXN TEMPLATE COMPLIANCE AND ENUMERATE CHEMICAL FRAGMENTS
polymerize = PolymerBuildProject.make_group(name='polymerize')

@polymerize
@PolymerBuildProject.pre(chemistry_valid)
@PolymerBuildProject.post(reactant_order_evaluated)
@PolymerBuildProject.operation
def determine_reactant_order(job : Job) -> None:
    '''
    Check that SMILES monomers are compatible with the 
    reaction template for the mechanism they claim to follow
    '''
    # 1) check that monomers fit a reaction template
    reactant_mol = load_job_rdmol(job)
    reactants = Chem.GetMolFrags(reactant_mol, asMols=True)
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
            logger.error(f'No valid ordering of reactants could be solved for the chosen "{job.doc.mechanism}" rxn template')

@polymerize
@PolymerBuildProject.pre(chemistry_valid) # TODO: find way to cache this from prior reactant order determination step
@PolymerBuildProject.pre(matches_rxn_template)
@PolymerBuildProject.post(functionalities_evaluated)
@PolymerBuildProject.operation
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
            reactant_mol = Chem.MolFromSmiles(reactant_smiles, sanitize=False) # CRITICAL that sanitize=False to avoid stripping
            Chem.SanitizeMol(reactant_mol, sanitizeOps=specification.SANITIZE_AS_KEKULE) # single, unified mol containing individual reactant as disconnected components
            
            num_funct_groups = substructures.num_substruct_queries_distinct(reactant_mol, rxn.GetReactantTemplate(i))
            functionalities.append(num_funct_groups)
            if num_funct_groups not in ALLOWED_FUNCTIONALITIES:
                logger.error(f'Found {num_funct_groups} active functional groups (vs any from {ALLOWED_FUNCTIONALITIES}) for molecule {reactant_smiles}')

        job.doc.functionalities = functionalities

@polymerize
@PolymerBuildProject.pre(matches_rxn_template)
@PolymerBuildProject.pre(monomers_satisfy_functionality)
@PolymerBuildProject.post(has_chemical_fragments)
@PolymerBuildProject.operation
def enum_fragments(job : Job) -> None:
    '''Enumerate all possible repeat unit fragment using cheminformatic reaction procedure'''
    rxn = load_job_rxn(job)
    reactor = PolymerizationReactor(rxn)

    reactant_mol = load_job_rdmol(job)
    reactants = Chem.GetMolFrags(reactant_mol, asMols=True)

    monogrp = MonomerGroup()
    with redirect_job_to_logfile(job) as logger:
        for intermediates, frags in reactor.propagate(reactants):
            for assoc_group_name, rdfragment in zip(ascii_uppercase, frags):
                # generate spec-compliant SMARTS
                raw_smiles = Chem.MolToSmiles(rdfragment)
                exp_smiles = specification.expanded_SMILES(raw_smiles)
                spec_smarts = specification.compliant_mol_SMARTS(exp_smiles)

                # record to monomer group
                affix = 'TERM' if MonomerGroup.is_terminal(rdfragment) else 'MID'
                monogrp.monomers[f'{assoc_group_name}_{affix}'] = [spec_smarts]

        monogrp.to_file(job.fn(PolymerBuildProject.FRAGMENTS_PATH))
        logger.info('Successfully enumerated and cached repeat unit fragments')

# 2) BUILD POLYMER STRUCTURE
oligomerize = PolymerBuildProject.make_group(name='oligomerize')

@oligomerize
# @PolymerBuildProject.pre.copy_from(determine_reactant_order)
# @PolymerBuildProject.pre.copy_from(enum_fragments)
@PolymerBuildProject.pre(has_chemical_fragments)
@PolymerBuildProject.post(coordinates_generated)
@PolymerBuildProject.operation
def build_oligomer_pdb(job : Job) -> None:
    '''Generate coordinates and build oligomer PDB file using mBuild'''
    seq = 'BA' # hard-coded for now, plan to make more flexible in the future
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
            DOP=2*(1 + (job.sp.DOP - 1)/len(seq)), # formula to convert target DOP (considering an AB pair as a repeat unit) to effective DOP in builder
            sequence=seq,
            energy_minimize=True, # TODO: add master config option for energy minimization at project level
        )
        mbmol_to_openmm_pdb(job.fn(PolymerBuildProject.OLIGOMER_PDB), polymer)
        logger.info('Successfully generated PDB structure file')

@oligomerize
# @PolymerBuildProject.pre.copy_from(build_oligomer_pdb)
@PolymerBuildProject.pre(coordinates_generated)
@PolymerBuildProject.post(chemical_info_assigned)
@PolymerBuildProject.operation
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

@oligomerize
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.post.true('r_eff') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('n_atoms_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('elem_counts_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.post.true('molar_mass_oligomer') # this ought to be fine, as these values should never be Falsy
@PolymerBuildProject.operation
def summarize_oligomer(job : Job) -> None:
    '''Compute simple summarizing info about an oligomer which streamline lattice packing,
    namely the number of atoms, the distribution of elements and the effect (max) radius'''
    offmol = load_job_oligomer_molecule(job)
    job.doc.r_eff = rdmol_effective_radius(offmol.to_rdkit()) # TODO: this shouldn't require an RDKit conversion, just access tot he conformer
    job.doc.n_atoms_oligomer = offmol.n_atoms
    job.doc.elem_counts_oligomer = elem_counts(offmol)
    job.doc.molar_mass_oligomer = sum(atom.mass for atom in offmol.atoms).magnitude
    
@oligomerize
@PolymerBuildProject.pre(chemical_info_assigned)
@PolymerBuildProject.post(partial_charges_assigned) # TODO: fill this in with something more substantive!!
@PolymerBuildProject.operation
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

# 3) PACK LATTICE
pack_lattice = PolymerBuildProject.make_group(name='pack_lattice') 

@pack_lattice
@PolymerBuildProject.pre(chemical_info_assigned) # don't need charges, only valid cornformer to pick sites
@PolymerBuildProject.post.true('n_oligomers')
@PolymerBuildProject.post.true('lattice_shape')
@PolymerBuildProject.post(lattice_sites_determined)
@PolymerBuildProject.operation
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

@pack_lattice
@PolymerBuildProject.pre(partial_charges_assigned)
@PolymerBuildProject.pre(lattice_sites_determined)
@PolymerBuildProject.post(neat_melt_packed)
@PolymerBuildProject.operation
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

@pack_lattice
@PolymerBuildProject.pre(neat_melt_packed)
@PolymerBuildProject.post(pbcs_determined)
@PolymerBuildProject.operation
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

# 4) PREPARE AND SERIALIZE OpenFF INTERCHANGE
to_interchange = PolymerBuildProject.make_group(name='to_interchange') 

@to_interchange
@PolymerBuildProject.pre(partial_charges_assigned)
@PolymerBuildProject.pre(neat_melt_packed)
@PolymerBuildProject.pre(pbcs_determined)
@PolymerBuildProject.pre(forcefield_is_valid)
@PolymerBuildProject.pre.not_(interchange_stereo_inconsistent)
@PolymerBuildProject.post(has_interchange)
@PolymerBuildProject.operation
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

# 5) EXPORT INTERCHANGE TO MD ENGINE FILES OF CHOICE
md_export = PolymerBuildProject.make_group(name='md_export') 

## LAMMPS versions
@PolymerBuildProject.label
def exported_to_LAMMPS(job : Job) -> bool:
    '''Check if LAMMPS files have been generated'''
    return all(
        has_nonempty_file(job, lmp_file)
            for lmp_file in PolymerBuildProject.LAMMPS_PATHS
    )

@md_export
@PolymerBuildProject.pre(has_interchange)
@PolymerBuildProject.post(exported_to_LAMMPS)
@PolymerBuildProject.operation
def export_LAMMPS_files(job : Job) -> None:
    '''Write LAMMPS data and input files for the melt topology'''
    interchange = load_job_interchange(job)
    Path(job.fn(PolymerBuildProject.LAMMPS_DIR)).mkdir(exist_ok=True) # ensure the child directory exists

    with redirect_job_to_logfile(job) as logger:
        logger.info('Writing LAMMPS input and data files')
        interchange_to_lammps(
            interchange=interchange,
            lmp_data_path=job.fn(PolymerBuildProject.LAMMPS_DATA_PATH),
            lmp_input_path=job.fn(PolymerBuildProject.LAMMPS_INPUT_PATH),
            lmp_data_filestr=f'"./{PolymerBuildProject.LAMMPS_DATA_PATH}"', # make this relative for portability
        )
        logger.info('LAMMPS files successfully written')

@PolymerBuildProject.pre(exported_to_LAMMPS)
@PolymerBuildProject.pre.never
@PolymerBuildProject.post.isfile(f'energies_LAMMPS.json')
@PolymerBuildProject.operation
def evaluate_energies_LAMMPS(job : Job) -> None:
    '''Evaluate starting structure energies of LAMMPS MD files'''
    ...

## OpenMM versions
@PolymerBuildProject.label
def exported_to_OpenMM(job : Job) -> bool:
    '''Check if OpenMM files have been generated'''
    return all(
        has_nonempty_file(job, lmp_file)
            for lmp_file in PolymerBuildProject.OPENMM_PATHS
    )

@md_export
@PolymerBuildProject.pre(has_interchange)
@PolymerBuildProject.pre.never
@PolymerBuildProject.post(exported_to_OpenMM)
@PolymerBuildProject.operation
def export_OpenMM_files(job : Job) -> None:
    '''Write OpenMM data and input files for the melt topology'''
    interchange = load_job_interchange(job)
    Path(job.fn(PolymerBuildProject.OPENMM_DIR)).mkdir(exist_ok=True) # ensure the child directory exists

    with redirect_job_to_logfile(job) as logger:
        integrator = VerletIntegrator(1*femtosecond) # hard-coded Intergrator now, will add support for more targetted integrator later
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

@PolymerBuildProject.pre(exported_to_LAMMPS)
@PolymerBuildProject.pre.never
@PolymerBuildProject.post.isfile(f'energies_OpenMM.json')
@PolymerBuildProject.operation
def evaluate_energies_OpenMM(job : Job) -> None:
    '''Evaluate starting structure energies of OpenMM MD files'''
    ...


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
        '-qp',
        '--quantity-precision',
        type=int,
        default=4,
        help='The number of decimal places to which to display and log physcal and numeric quantities',
    )
    parser.add_argument(
        '--strict-stereo',
        action='store_true',
        help='Optional, whether to enforce strict and unambiguous stereochemistry when loading molecules into OpenFF toolkit objects',
    )
    # TODO : implement log level setting

    # separate this scripts args from those required by signac
    start_args, signac_args = parser.parse_known_args()
    assert start_args.project_path.exists() and start_args.project_path.is_dir()
    
    # configure global vars in Project definition and initialize project instance
    PolymerBuildProject.QUANTITY_PRECISION = start_args.quantity_precision
    PolymerBuildProject.RELAXED_STEREO = not start_args.strict_stereo
    new_project = PolymerBuildProject.get_project(start_args.project_path)

    # mock remaining Signac args for parser and run project CLI interface
    sys.argv[1:] = signac_args # NOTE: this is an ugly hack to allow this script to take CLI args while not disturbing Signacs tastes for arguments
    new_project.main()

if __name__ == '__main__':
    main()
