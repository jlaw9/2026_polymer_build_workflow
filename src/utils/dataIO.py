'''For reading and writing monomer data files'''

import logging
from typing import Container, Iterable, Optional

import json
import pandas as pd
from pathlib import Path


# VALIDATNG FILES
def validate_file_path(
        path : Path,
        check_missing : bool=False,
        check_already_exists : bool=False,
        check_has_extension : bool=False,
        valid_extensions : Optional[Container[str]]=None,
    ) -> None:
    '''Check that a path exists, is pathlike, and has a valid extension
    Performs no check by default if no arguments other than the path are passed; these NEED to be supplied when called!'''
    # Meta-error for nonsensical argument checks
    if check_missing and check_already_exists:
        raise ValueError('Incongruent checks requested; invalid to check both that a file already exists AND is missing')
    
    # Conditional file checks
    if check_missing and not path.exists():
        raise FileNotFoundError(f'No file exists at "{path}"')
    if check_already_exists and path.exists():
        raise PermissionError(f'File already exists at "{path}"')
    if check_has_extension and not path.suffix:
        raise IsADirectoryError(f'Input file missing file extension, appears like directory ("{path}")')
    if valid_extensions and (path.suffix not in valid_extensions):
        raise ValueError(
            f'Cannot read data from {path.suffix} file, choose one of' \
            f'the following valid extensions: {[ext for ext in valid_extensions]}' # NOTE: using comprehension instead of list() to give expected output for dicts
        )

# READING DATA FOR PROJECT SETUP
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