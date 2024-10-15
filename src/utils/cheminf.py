'''Cheminformatic and SMILES/SMARTS string-related functionality'''

import logging

from typing import Optional, Sequence, Union
from ast import literal_eval

from rdkit.Chem import CanonSmiles
from polymerist.smileslib.primitives import is_valid_SMILES


def parse_monomer_smiles(smiles : Union[str, Sequence[str]], canonicalize : bool=True) -> Optional[str]:
    '''Enforces formatting of SMILES monomer inputs as a single dot-bond joined string of monomer-wise SMILES'''
    if isinstance(smiles, str): # convert Sequences saved as strings to literal Sequences (i.e.  "('A', 'B')" -> ('A', 'B'))
        try:
            smiles = literal_eval(smiles)
        except (SyntaxError, ValueError):
            logging.debug(f'SMILES stayed as {smiles} (no tuple-ification detected)')
    
    # print('#', smiles)
    if isinstance(smiles, Sequence) and not isinstance(smiles, str): # strings are technically Sequences, but we don't want to reformat them here
        smiles = '.'.join(smiles)

    # print('##', smiles)
    if not (isinstance(smiles, str) and is_valid_SMILES(smiles)):
        # raise TypeError
        return None
    
    if canonicalize:
        logging.info(f'Canonicalizing SMILES string "{smiles}"')
        smiles = CanonSmiles(smiles)
    
    return smiles