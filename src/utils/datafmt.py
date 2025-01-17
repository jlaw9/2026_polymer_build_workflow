'''For formatting and stripping fields out of monomer databases'''

import logging

from typing import Iterable, TypeAlias
StringMap : TypeAlias = dict[str, str]

import re
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

def standardize_monomer_data_columns(dataframe : pd.DataFrame) -> None:
    '''Standardize column naming and format of required monomer data DataFrame (in-place)'''
    STATEPOINT_ATTR_COLUMNS : dict[str, tuple[str]] = { # the attributes to save and the column(s) to check for these values
        'smiles_original' : ('smiles_monomer', 'monomer', 'monomers', 'Monomer', 'Monomers'),
        'mechanism' : ('mechanism', 'rxnname', 'Chemistry')
    }
    attr_locs = locate_attr_cols(dataframe, STATEPOINT_ATTR_COLUMNS) # this will raise Exception if any of the fields cannot be found
    dataframe.rename(
        columns={col_found_in : std_name for std_name, col_found_in in attr_locs.items()},
        inplace=True # perform rename in-place to avoid allocating memory for new (potentially large) dataframe
    )

def label_monomer_statepoint_data(dataframe : pd.DataFrame, rxn_mapping : dict[str, str], uniquify_chemistry : bool=False) -> None:
    '''Ensures monomer data DataFrame has detailed SMILES and reaction info columns,
    and that metadata vs and essential statepoint data column labels are clearly differentiated
    '''
    logging.info('Canonicalizing all SMILES')
    dataframe['smiles_canonical'] = dataframe['smiles_original'].map(lambda smi : parse_monomer_smiles(smi, canonicalize=True))
    if uniquify_chemistry:
        logging.info('Purging monomer records with duplicate chemistries')
        dataframe.drop_duplicates('smiles_canonical', inplace=True)
    
    logging.info('Expanding SMILES to be chemically explicit')
    dataframe['smiles_explicit' ] = dataframe['smiles_canonical'].map(lambda smi : expanded_SMILES(smi, assign_map_nums=False, kekulize=False))
    
    logging.info('Looking up reaction mechanism')
    dataframe['rxn_smarts'] = dataframe['mechanism'].map(rxn_mapping)

    statepoint_colnames = ['smiles_explicit', 'rxn_smarts'] # explicitly mark statepoint and metadata columns to simplify final parse
    dataframe.rename(
        columns=lambda colname : f'{colname}<statedata>' if colname in statepoint_colnames else f'{colname}<metadata>',
        inplace=True # perform rename in-place to avoid allocating memory for new (potentially large) dataframe 
    ) 

def parse_field_names_and_roles(dataframe : pd.DataFrame) -> tuple[StringMap, StringMap]:
    '''Extract the field (column) names and data role metadata
    Returns two dicts mapping from column names as-they-are to names and data roles, respectively'''

    HEADER_ROLE_RE = re.compile('(?P<field_name>.*?)<(?P<field_role>.*?)>') # role is delimited by chevrons

    colname_tag_free  : StringMap = {}
    colname_data_role : StringMap = {}
    for colname in dataframe.columns:
        matches = re.match(HEADER_ROLE_RE, colname)
        assert matches is not None

        match_fields : StringMap = matches.groupdict()
        colname_tag_free[colname] = match_fields['field_name']
        colname_data_role[colname] = match_fields['field_role']
        
    return colname_tag_free, colname_data_role