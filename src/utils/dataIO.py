'''For reading and writing monomer data files'''

import logging
from typing import Iterable

import json
import pandas as pd
from pathlib import Path

from .filelib import validate_file_path


READER_FNS_BY_EXT = {
    '.xlsx' : pd.read_excel,
    '.csv'  : pd.read_csv,
}
WRITER_FNS_BY_EXT = {
    '.xlsx' : pd.DataFrame.to_excel,
    '.csv'  : pd.DataFrame.to_csv,
}
RXN_MAP_EXTS = ('.json',)

def read_monomer_data(mdat_paths : Iterable[Path]) -> list[pd.DataFrame]:
    '''Validate and read in a series of monomer data paths
    Returns a list of dataframes, containing monomer data in the order that paths were passed'''
    mdat_dataframes : list[pd.DataFrame] = []
    for mdat_path in mdat_paths:
        validate_file_path(mdat_path, valid_extensions=READER_FNS_BY_EXT, check_missing=True)
        reader_fn = READER_FNS_BY_EXT[mdat_path.suffix] # don't use get() here; WANT a KeyError if invalid
        logging.info(f'Reading monomer data from {mdat_path}')
        mdat_dataframes.append(reader_fn(mdat_path))

    return mdat_dataframes

def read_rxn_mapping_data(rxn_mapping_path : Path) -> dict[str, str]:
    '''Validate and read in a reaction mapping datafile
    Returns a dict keyed by reaction anem whose keys are SMARTS for the corresponding functional reaction'''
    validate_file_path(rxn_mapping_path, valid_extensions=RXN_MAP_EXTS, check_missing=True)
    with rxn_mapping_path.open('r') as rxn_map_file:
        logging.info(f'Reading reaction data from {rxn_mapping_path}')
        return json.load(rxn_map_file)