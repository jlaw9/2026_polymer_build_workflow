'''For formatting and stripping fields out of monomer databases'''

import logging

from typing import Iterable
import pandas as pd

from .containers import stringify_dict
from .cheminf import parse_monomer_smiles

from polymerist.polymers.monomers.specification import expanded_SMILES


def locate_attr_cols(dataframe : pd.DataFrame, columns_to_check : dict[str, Iterable[str]]) -> dict[str, str]:
    '''Takes a dataframe of monomer training data and a dict of desired attributes and the columns in the dataframe it might be found in
    Checks that those columns are present and returns dict with first column for each if all are present, or NoneType otherwise'''
    attr_columns : dict[str, str] = {}
    for targ_attr, col_names_to_check in columns_to_check.items():
        for col_name in col_names_to_check:
            if col_name in dataframe:
                attr_columns[targ_attr] = col_name
                break
        else:
            raise IndexError(f'No matching columns for attribute "{targ_attr} were found from queries: "{col_names_to_check}"')
    logging.info('Found valid columns name mappings:\n\t' + stringify_dict(attr_columns))
        
    return attr_columns

def standardize_monomer_data(dataframe : pd.DataFrame, rxn_mapping : dict[str, str]) -> None:
    '''
    Standardize column naming and format of required monomer data DataFrame (in-place)
    and ensure required fields for statepoint generation are present
    '''
    STATEPOINT_ATTR_COLUMNS : dict[str, tuple[str]] = { # the attributes to save and the column(s) to check for these values
        'smiles_original' : ('smiles_monomer', 'monomer', 'monomers', 'Monomer', 'Monomers'),
        'mechanism' : ('mechanism', 'rxnname', 'Chemistry')
    }
    attr_locs = locate_attr_cols(dataframe, STATEPOINT_ATTR_COLUMNS) # this will raise Exception if any of the fields cannot be found
    dataframe.rename(
        columns={col_found_in : std_name for std_name, col_found_in in attr_locs.items()},
        inplace=True # perform rename in-place to avoid allocating memory for new (potentially large) dataframe
    )

    # insert new columns for processed statepoint data
    logging.info('Canonicalizing all SMILES')
    dataframe['smiles_canonical'] = dataframe['smiles_original'].map(lambda smi : parse_monomer_smiles(smi, canonicalize=True))
    
    logging.info('Expanding SMILES to be chemically explicit')
    dataframe['smiles_explicit' ] = dataframe['smiles_canonical'].map(lambda smi : expanded_SMILES(smi, assign_map_nums=False))
    
    logging.info('Looking up reaction mechanism')
    dataframe['rxn_smarts'] = dataframe['mechanism'].map(rxn_mapping)

    statepoint_colnames = ['smiles_explicit', 'rxn_smarts'] # explicitly mark statepoint and metadata columns to simplify final parse
    dataframe.rename(
        columns=lambda colname : f'{colname}<statedata>' if colname in statepoint_colnames else f'{colname}<metadata>',
        inplace=True # perform rename in-place to avoid allocating memory for new (potentially large) dataframe 
    ) 

    return dataframe