'''For cleaning up and standardizing raw monomer input data files shipped from NREL'''

# Logging
import logging

import warnings
warnings.filterwarnings(action='ignore')

# Stdlib
from typing import Any, Container, Iterable, Sequence, Optional, Union

from ast import literal_eval
from argparse import ArgumentParser, Namespace

# File I/O
import json
from pathlib import Path

import pandas as pd

# Cheminformatics
from rdkit.Chem import CanonSmiles

# Custom
from polymerist.smileslib.primitives import is_valid_SMILES
from polymerist.polymers.monomers import specification


# VALIDATING PATHS
def _stringify_dict(anydict : dict[Any, Any], sep : str=', ', joiner : str='-->') -> str:
    '''Create an inline string describing the mappping in a dict'''
    return sep.join(
        f'({key!s} {joiner} {value!s})'            
            for key, value in anydict.items()
    )

def _validate_file_path(
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
    
# SMILES PARSING
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
    
# DATAFRAME MANIPULATION
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
    logging.info('Found valid columns name mappings:\n\t' + _stringify_dict(attr_columns))
        
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
    dataframe['smiles_explicit' ] = dataframe['smiles_canonical'].map(lambda smi : specification.expanded_SMILES(smi, assign_map_nums=False))
    
    logging.info('Looking up reaction mechanism')
    dataframe['rxn_smarts'] = dataframe['mechanism'].map(rxn_mapping)

    statepoint_colnames = ['smiles_explicit', 'rxn_smarts'] # explicitly mark statepoint and metadata columns to simplify final parse
    dataframe.rename(
        columns=lambda colname : f'{colname}<S>' if colname in statepoint_colnames else f'{colname}<M>',
        inplace=True # perform rename in-place to avoid allocating memory for new (potentially large) dataframe 
    ) 

    return dataframe

# READING INPUT DATA
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
        _validate_file_path(mdat_path, valid_extensions=READER_FNS_BY_EXT, check_missing=True)
        reader_fn = READER_FNS_BY_EXT[mdat_path.suffix] # don't use get() here; WANT a KeyError if invalid
        logging.info(f'Reading monomer data from {mdat_path}')
        mdat_dataframes.append(reader_fn(mdat_path))

    return mdat_dataframes

def read_rxn_mapping_data(rxn_mapping_path : Path) -> dict[str, str]:
    '''Validate and read in a reaction mapping datafile
    Returns a dict keyed by reaction anem whose keys are SMARTS for the corresponding functional reaction'''
    _validate_file_path(rxn_mapping_path, valid_extensions=RXN_MAP_EXTS, check_missing=True)
    with rxn_mapping_path.open('r') as rxn_map_file:
        logging.info(f'Reading reaction data from {rxn_mapping_path}')
        return json.load(rxn_map_file)
    
def sanitize_monomer_data_paths(args : Namespace) -> list[Path]:
    '''Handles both the direct "monomer-paths" or "glob" modes of passing monomer input files'''
    if not args.output_dir.is_dir():
        args.output_dir.mkdir()

    if args.glob is None:
        return args.monomer_paths
    elif args.monomer_paths is None:
        return [path for path in Path.cwd().glob(args.glob)]

    
# OPERATION MODES
def format_merged(args : Namespace) -> None:
    '''
    Takes one or more monomer data files, reformats them, merges them into a single
    "master" file with shared columns, and writes this to a single output file
    '''
    output_path = args.output_dir / args.output_file
    _validate_file_path(output_path, check_missing=False, check_already_exists=not args.allow_overwrites, valid_extensions=WRITER_FNS_BY_EXT)

    rxn_mapping = read_rxn_mapping_data(args.rxn_mapping_path)
    monomer_paths = sanitize_monomer_data_paths(args)
    monomer_dfs = read_monomer_data(monomer_paths)
    for df in monomer_dfs:
        standardize_monomer_data(df, rxn_mapping=rxn_mapping)

    master_df = pd.concat(monomer_dfs) # columns we care about should be aligned now that the dataframes are standardized
    if args.keep_n is not None:
        keep_n = min(args.keep_n, len(master_df)) # clamp number of sample to the size of the dataset
        if args.random:
            logging.info(f'Sampling a random {keep_n} records from merged dataset')
            master_df = master_df.sample(keep_n) # randomly sample N records
        else:
            logging.info(f'Sampling the first {keep_n} records from merged dataset')
            master_df = master_df.head(keep_n) # sample the first N records, unless N exceed the size of the database
    
    writer_fn = WRITER_FNS_BY_EXT[output_path.suffix]
    logging.info(f'Writing monomer data to {output_path}...')
    writer_fn(master_df, output_path) # need to call with master_df in place of "self" argument
    logging.info(f'Write to {output_path} complete')

def format_sequential(args : Namespace) -> None:
    '''
    Takes one or more monomer data files, individually reformats them and
    writes the reformatted files to the file destination(s) of choice

    In either the "postfix" or "new-names" modes, each outfile file will 
    have the same extension as its respective input file
    '''
    monomer_paths = sanitize_monomer_data_paths(args)
    
    # determine output file paths - argparse check ensures the postfix and new_names modes will have mutually exclusive options
    if args.new_names is None: 
        logging.info('Determining output file names via the "postfix" directive')
        output_paths : list[Path] = []
        for input_path in monomer_paths:
            output_path = args.output_dir / f'{input_path.stem}_{args.postfix}{input_path.suffix}'
            _validate_file_path(output_path, check_missing=False, check_already_exists=not args.allow_overwrites, valid_extensions=WRITER_FNS_BY_EXT)
            output_paths.append(output_path)

    elif args.postfix is None: 
        logging.info('Determining output file names via the "new_names" directive')
        if (num_names := len(args.new_names)) != (num_inputs := len(monomer_paths)):
            raise IndexError(
                'If providing output names, must provide exactly as many names as monomer input files\n' \
                f'Provides {num_inputs} data files, but {num_names} output names'
            )
        
        output_paths : list[Path] = []
        for output_name, input_path in zip(args.new_names, monomer_paths):
            output_path = args.output_dir / f'{output_name}{input_path.suffix}'
            _validate_file_path(output_path, check_missing=False, check_already_exists=not args.allow_overwrites, valid_extensions=WRITER_FNS_BY_EXT)
            output_paths.append(output_path)

    # read, reformat, and write out data
    rxn_mapping = read_rxn_mapping_data(args.rxn_mapping_path)
    monomer_dfs = read_monomer_data(monomer_paths)

    for monomer_df, output_path in zip(monomer_dfs, output_paths):
        standardize_monomer_data(monomer_df, rxn_mapping=rxn_mapping)
        # NOTE: don't need to validate output, as this was done in the output path compile step prior
        writer_fn = WRITER_FNS_BY_EXT[output_path.suffix]
        logging.info(f'Writing monomer data to {output_path}...')
        writer_fn(monomer_df, output_path) # need to call with master_df in place of "self" argument
        logging.info(f'Write to {output_path} complete')


# ARGUMENT PARSING AND DISPATCH
def main() -> None:
    '''Read monomer input data, standardize column data, and ensure required fields for statepoint generation are present'''
    parser = ArgumentParser(description='Clean up and standardized NREL PolyID monomer training data files for MD structure build')
    subparsers = parser.add_subparsers()

    # auxiliary parser for handling shared input parameters between both substrategies
    input_parser = ArgumentParser(add_help=False) 
    file_input_group = input_parser.add_mutually_exclusive_group(required=True)
    file_input_group.add_argument(
        '-mdat',
        '--monomer-paths',
        type=Path,
        nargs='+',
        help='The relative path(s) to the tabular data file(s) containing monomer SMILES and other training data',
    )
    file_input_group.add_argument(
        '-g',
        '--glob',
        help='A glob pattern to use to locate input files',
    )

    input_parser.add_argument( 
        '-rxns',
        '--rxn_mapping-path',
        type=Path,
        required=True,
        help='The path to a JSON file containing a reaction name mapping\n' \
            'Should contain dict whose keys are reaction names and whose values are reaction SMARTS strings'
    )
    input_parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='The directory into which the individual formatted datafiles should be saved to'
    )
    input_parser.add_argument(
        '-aow',
        '--allow-overwrites',
        action='store_true',
        help='Whether to permit overwriting output files which already exist (default is False)'
    )

    # subparsers for different modes of operation
    parser_merge = subparsers.add_parser('merge', parents=[input_parser], description=format_merged.__doc__)
    parser_merge.add_argument(
        '-of',
        '--output-file',
        type=Path,
        required=True,
        help='File path that the single, merged and formatted dataset should be written to'
    )
    parser_merge.add_argument(
        '-k',
        '--keep-n',
        type=int,
        help='Optional number of records to subsample from the dataset, if the full dataset is not desired'
    )
    parser_merge.add_argument(
        '-rand',
        '--random',
        action='store_true',
        help='When a value is set for --keep-n, dictates whether the subsample should be the first keep_n (default) or a random sample of keep_n'
    )
    parser_merge.set_defaults(func=format_merged)

    parser_sequence = subparsers.add_parser('sequence', parents=[input_parser], description=format_sequential.__doc__)
    rename_group = parser_sequence.add_mutually_exclusive_group(required=True)
    rename_group.add_argument(
        '-pf',
        '--postfix',
        help='An optional postfix to append to the name of the original file when saving processed output'
    )
    rename_group.add_argument(
        '-nn',
        '--new-names',
        nargs='+',
        help='The new name(s) (including extension!) to assign to the output file\n' \
            'Must NOT match any of the input filenames!'
    )
    parser_sequence.set_defaults(func=format_sequential)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO) # TODO: change log level with args
    # logging.basicConfig(level=logging.ERROR) # TODO: change log level with args
    args.func(args)

if __name__ == '__main__':
    main()