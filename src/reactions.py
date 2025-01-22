'''Define and cache reaction templates for polymerization procedure'''

import logging
logging.basicConfig(level=logging.INFO)

from typing import Sequence, Optional

from pathlib import Path

import json
from dataclasses import dataclass, field

from rdkit import Chem

from polymerist.genutils.fileutils.jsonio.jsonify import make_jsonifiable
from polymerist.polymers.monomers.specification import expanded_SMILES

from polymerist.rdutils.bonding import portlib
from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.reactions.assembly import ReactionAssembler

try: # call as python module
    from .__init__ import _parent_dir
    from .utils.logs import format_error_for_log
except ImportError: # call as script file
    from __init__ import _parent_dir
    from utils.logs import format_error_for_log


# INITIALIZING SMARTS FOR FUNCTIONAL GROUPS
fn_group_smarts = { # mapped SMARTS (SMIRKS) for common functional groups
    'phthalimide'     : '[*:1]-[C:2](=[O:3])-[N:4](-[*:5])-[C:6](=[O:7])-[*:8]',
    'carbamate'       : '[*:1]-[N:2](-[C:3](=[O:4])-[O:5]-[*:6])-[H:7]',
    'cyclocarbonate'  : '[*:1]-[C:2]1(-[H:8])-[O:3]-[C:4](=[O:5])-[O:6]-[C:7]-1(-[H:9])-[H:10]',
    'anhydride'       : '[*:1]-[C:2](=[O:3])-[O:4]-[C:5](=[O:6])-[*:7]',
    'vinyl'           : '[*:1]-[C:2](=[C:3](-[H:5])-[H:6])-[H:4]',
    'terminal_alkene' : '[*:1]-[C:2](-[*:3])=[C:4](-[H:5])-[H:6]',
    'acyl_chloride'   : '[Cl:1]-[C:2](=[O:3])-[*:4]',
    'carboxyl'        : '[O:1](-[C:2](=[O:3])-[*:4])-[H:5]',
    'ester'           : '[*:1]-[O:2]-[C:3](=[O:4])-[*:5]',
    'amine'           : '[N:1](-[*:2])(-[H:3])-[H:4]',
    'hydroxyl'        : '[O:1](-[*:2])-[H:3]',
    # 'hydroxyl'        : '[O:1](-[C:2]-[*:3])-[H:4]',
    'isocyanate'      : '[O:1]=[C:2]=[N:3]-[*:4]'
}
fn_group_mols : dict[str, Chem.Mol] = {}
for group_name, smarts in fn_group_smarts.items():
    fn_group_mol = Chem.MolFromSmiles(smarts, sanitize=False) # despite being SMARTS, the molecule needs to be initialized via SMILES to avoid query persistence during rearrangements
    fn_group_mol.SetProp('_Name', group_name) # shows up as MDL mol label in file
    
    group_label_assigned : bool = False # ensure this is only done for one linker per group
    for linker_atom_id in portlib.get_linker_ids(fn_group_mol):
        linker_atom = fn_group_mol.GetAtomWithIdx(linker_atom_id)
        linker_atom.SetProp('dummyLabel', 'R') # label linkers as R-groups
        if not group_label_assigned:
            linker_atom.SetProp('molFileValue', group_name)
            group_label_assigned = True
    
    fn_group_mols[group_name] = fn_group_mol

# INITIALIZING BYPRODUCTS TEMPLATES
byproduct_smarts = { # unmapped SMARTS for extraneous molecules which are to be unmapped and removed
    'water'   : '[H]-[O]-[H]',
    'hcl'     : '[H]-[Cl]',
    'alcohol' : '[*]-[O]-[H]'
}
byproduct_mols = {
    byprod_name : Chem.MolFromSmarts(smarts)
        for byprod_name, smarts in byproduct_smarts.items()
}

# DEFINING MAPPING FROM POLYID MECHANISM NAMES TO MORE DESCRIPTIVE NAMES USED HERE
polyid_backmap = {
    "amide"     : "polyamide",
    "carbonate" : "polycarbonate_phosgene",
    "ester"     : "polyester",
    "imide"     : "polyimide",
    "urethane"  : "polyurethane_isocyanate",
    "NIPU"      : "polyurethane_nonisocyanate",
    "vinyl"     : "polyvinyl_head_tail"
}

# DEFINING PARAMETERS FOR REACTION TEMPLATES 
@make_jsonifiable
@dataclass
class ReactionInfo:
    '''Bundles together minimal information needed to assemble a reaction template'''
    reactant_groups     : Sequence[str]
    byproduct_templates : Sequence[str]
    bond_derangement : dict[int, tuple[int, int]]
    test_reactant_smiles : Sequence[str]

rxn_inputs : dict[str, ReactionInfo] = {
    'polyester' : ReactionInfo(
        reactant_groups=['hydroxyl', 'carboxyl'],
        byproduct_templates=['water'],
        bond_derangement={
            1 : (3, 5),
            4 : (5, 3)
        },
        test_reactant_smiles=('OCCO', 'O(C=O)c1ccc(cc1)C(=O)O'), # PET,
    ),
    'polyamide' : ReactionInfo(
        reactant_groups=['amine', 'carboxyl'],
        byproduct_templates=['water'],
        bond_derangement={
            1 : (3, 6),
            5 : (6, 3)
        },
        test_reactant_smiles=('NCCCCCCN', 'O=C(O)CCCCC(=O)O'), # Nylon-6,6
    ),
    'polyimide' : ReactionInfo(
        reactant_groups=['amine', 'anhydride'],
        byproduct_templates=['water'],
        bond_derangement={
            4 : (1, 8),
            6 : (8, 1),
            3 : (1, 8), # doubles up carbonyl transfer - must have target atoms as beginning to maintain canonical derangement form
            9 : (8, 1),
        },
        test_reactant_smiles=('O(c1ccc(N)cc1)c2ccc(cc2)N', 'C1=C2C(=CC3=C1C(=O)OC3=O)C(=O)OC2=O'), # DuPont Kapton (poly (4,4'-oxydiphenylene-pyromellitimide))
    ),
    'polycarbonate_phosgene' : ReactionInfo(
        reactant_groups=['hydroxyl', 'acyl_chloride'],
        byproduct_templates=['hcl'],
        bond_derangement= {
            1 : (3, 5),
            4 : (5, 3)
        },
        test_reactant_smiles=('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'ClC(=O)Cl'), # BPA + phosgene
    ),
    'polycarbonate_nonphosgene' : ReactionInfo(
        reactant_groups=['hydroxyl', 'ester'],
        byproduct_templates=['alcohol'],
        bond_derangement= {
            5 : (6, 3),
            1 : (3, 6)
        },
        test_reactant_smiles=('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'O=C(Oc1ccccc1)Oc2ccccc2'), # BPA + diphenyl carbonate
    ),
    'polyurethane_isocyanate' : ReactionInfo(
        reactant_groups=['isocyanate', 'hydroxyl'],
        byproduct_templates=[],
        bond_derangement={
            7 : (5, 3),
            2 : (3, 5)
        },
        test_reactant_smiles=('CC(=C)C(=O)OCC1COC(=O)O1', 'NCCCCCCN'), # PCA (propylene carbonate acrylate) + hexamethylenediamine
    ),
    'polyurethane_nonisocyanate' : ReactionInfo(
        reactant_groups=['cyclocarbonate', 'amine'],
        byproduct_templates=[],
        bond_derangement= {
            5  : (4, 11), # (7, 11)
            13 : (11, 4)  # (11, 7)
        },
        test_reactant_smiles=('O=C=N\CCCCCC/N=C=O', 'OCCCCO'), # Bayer HDI + BDO
    ),
    'polyvinyl_head_tail' : ReactionInfo(
        reactant_groups=['terminal_alkene', 'vinyl'],
        byproduct_templates=[],
        bond_derangement= {
            4  : (2, 8),
            12 : (8, 2)
        },
        test_reactant_smiles=('c1ccccc1C=C', 'c1ccccc1C=C'), # polystyrene
    ),
}

# INITIALIZING REACTION ASSEMBLERS (WITHOUT INITIATING ASSEMBLY PROCESS)
rxn_assemblers : dict[str, ReactionAssembler] = {}
test_reactants_catalogue : dict[str, list[Chem.Mol]] = {}

num_rxn_inputs = len(rxn_inputs)
for i, (rxnname, rxninfo) in enumerate(rxn_inputs.items(), start=1):
    logging.info(f'Initializing reaction template {i}/{num_rxn_inputs} ("{rxnname}")')
    rxn_assembler = ReactionAssembler(
        reactive_groups=[fn_group_mols[reacgrp_name] for reacgrp_name in rxninfo.reactant_groups],
        byproducts=[byproduct_mols[byprod_name] for byprod_name in rxninfo.byproduct_templates],
        bond_derangement=rxninfo.bond_derangement,
    )
    rxn_assemblers[rxnname] = rxn_assembler # this is solely for debug in external modules

    logging.info('Initializing test reactants for validation')
    test_reactants = []
    for smiles in rxninfo.test_reactant_smiles:
        exp_smiles = expanded_SMILES(smiles, assign_map_nums=False, kekulize=False)
        reactant_mol = Chem.MolFromSmiles(exp_smiles, sanitize=False)
        Chem.SanitizeMol(reactant_mol) # implicitly undoes kekulization anyway
        test_reactants.append(reactant_mol)
    test_reactants_catalogue[rxnname] = [Chem.Mol(reactant) for reactant in test_reactants]
    

if __name__ == '__main__':
    RXNS_DIR = _parent_dir / 'reactions'
    RXNS_DIR.mkdir(exist_ok=True)

    # ASSEMBLING AND TESTING REACTIONS
    rxns : dict[str, AnnotatedReaction] = {}
    rxn_smarts : dict[str, str] = {}

    num_rxn_inputs = len(rxn_inputs)
    for i, (rxnname, rxninfo) in enumerate(rxn_inputs.items(), start=1):
        logging.info(f'Assembling reaction {i}/{num_rxn_inputs} ("{rxnname}")')
        rxn_inputs_path = RXNS_DIR / f'{rxnname}_inputs.json'
        rxninfo.to_file(rxn_inputs_path)

        logging.info('Initializing Reaction')
        rxn_assembler = rxn_assemblers[rxnname]
        rxn = rxn_assembler.assemble_rxn(show_steps=False)
        rxn.rxnname = rxnname
        logging.info('Reaction successfully assembled')

        logging.info('Validating reaction template on test reactants')
        test_reactants = test_reactants_catalogue[rxnname]
        try:
            reactor = PolymerizationReactor(rxn)
            _ = [(dimer, frags) for dimer, frags in reactor.propagate(test_reactants)]
            logging.info('VALIDATION SUCCESSFUL: Reaction and test reactants are compatible')
        except Exception as error: # TODO: make this more granular
            logging.error(format_error_for_log(error))
            continue

        logging.info('Recording reaction object and representative SMARTS')
        rxns[rxnname] = rxn
        rxn_smarts[rxnname] = rxn.to_smarts().replace('#0', '*')  # temporary fix to 0-atomic number bug
    
    # BUILDING DIRECT MAP FROM POLYID MECHANISM NAMES TO SMARTS FOR BREVITY
    rxns_polyid = { # expand this here to avoid need for backmap awareness downstream
        polyid_name : rxn_smarts[descriptive_name]
            for polyid_name, descriptive_name in polyid_backmap.items() 
    }
    
    # PATHS WHERE REACTION OUTPUTS SHOULD LIVE (HARD-CODED FOR NOW)
    rxn_info_paths : dict[Path, dict[str, str]] = {
        Path(RXNS_DIR, 'functional_groups.json') : fn_group_smarts,
        Path(RXNS_DIR, 'byproducts.json')        : byproduct_smarts,
        Path(RXNS_DIR, 'rxn_smarts.json')        : rxn_smarts,
        Path(RXNS_DIR, 'rxns_polyID.json')       : rxns_polyid,
        Path(RXNS_DIR, 'rxn_names_polyID.json')  : polyid_backmap,
        # Path(RXNS_DIR, 'test_reactants.json')    : test_reactants_catalogue, # NOTE: not serialized, since this dict contains Mol objects (not SMILES strings)
    }
    for target_path, rxn_info_dict in rxn_info_paths.items():
        # if not target_path.exists():
        with open(target_path, 'w') as file:
            json.dump(rxn_info_dict, file, indent=4)

    # for chemistry, rxn in rxns.items():
        # rxn.to_rxnfile(RXNS_DIR / f'{chemistry}.rxn', wilds_to_R_groups=True)