'''Cheminformatic and SMILES/SMARTS string-related functionality'''

import logging
from typing import Iterable, Optional, Sequence, Union

from pathlib import Path
from ast import literal_eval

from rdkit import Chem

from polymerist.smileslib.cleanup import is_valid_SMILES, expanded_SMILES
from polymerist.genutils.decorators.functional import allow_string_paths


@allow_string_paths
def is_valid_sdfile(
    sdfile_path : Path,
    sanitize : bool=False,
    remove_hs : bool=False,
    strict_parsing : bool=False,
) -> bool:
    '''Check if a chemical SDF contains all valid serialized molecules'''
    try:
        supplier = Chem.SDMolSupplier(
            sdfile_path,
            sanitize=sanitize,
            removeHs=remove_hs,
            strictParsing=strict_parsing
        )
    except OSError:
        return False

    try:
        for mol in supplier:
            if mol is None:
                return False
        else:
            return True
    except StopIteration:
        return False
    
def parse_monomer_smiles(smiles : Union[str, Sequence[str]], canonicalize : bool=True) -> Optional[str]:
    '''Enforces formatting of SMILES monomer inputs as a single dot-bond joined string of monomer-wise SMILES'''
    if isinstance(smiles, str): # convert Sequences saved as strings to literal Sequences (i.e.  "('A', 'B')" -> ('A', 'B'))
        try:
            smiles = literal_eval(smiles)
        except (SyntaxError, ValueError):
            logging.debug(f'SMILES stayed as {smiles} (no tuple-ification detected)')
    
    if isinstance(smiles, Sequence) and not isinstance(smiles, str): # strings are technically Sequences, but we don't want to reformat them here
        smiles = '.'.join(smiles)

    if not (isinstance(smiles, str) and is_valid_SMILES(smiles)):
        # raise TypeError
        return None
    
    if canonicalize:
        logging.info(f'Canonicalizing SMILES string "{smiles}"')
        smiles = Chem.CanonSmiles(smiles)
    
    return smiles
