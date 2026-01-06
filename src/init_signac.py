'''Initialize project workspace directory from statepoints defined by build project config and chemical input dataset'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
from argparse import ArgumentParser, Namespace

from typing import Optional

import json
from pathlib import Path
from shutil import copytree

import numpy as np
import pandas as pd

from signac import init_project

from polymerist.genutils.iteration import cartesian_grid

from . import _parent_dir
from .parameters import SystemParameters, PARAMS_SWEPT_PATH, standardize_params_swept
from .utils.dataIO import validate_file_path, read_monomer_data


# Helper functions
def read_monomer_dataframe(
    mono_data_path : Path,
    number_to_sample : Optional[int]=None,
    random : bool=False
) -> pd.DataFrame:
    '''Read and subselect monomer data fields according to provided read parameters'''
    monomer_df = read_monomer_data([mono_data_path])[0] # take first of one as THE dataframe
    monomer_df.set_index(monomer_df.columns[0], inplace=True) # alternative to index_cols arg in read function
    monomer_df.replace(np.nan, None, inplace=True) # convert NaN values to JSON-serializable NoneType
    
    if number_to_sample is not None:
        if random:
            monomer_df = monomer_df.sample(number_to_sample)
        else:
            monomer_df = monomer_df.head(min(number_to_sample, len(monomer_df)))

    return monomer_df

def generate_statepoints(args : Namespace) -> None:
    '''Initialize project directory and job statepoint files for chosen monomer data fields'''
    # read monomer data into memory
    monomer_df = read_monomer_dataframe(
        args.monomer_data,
        number_to_sample=args.number_to_sample,
        random=args.random
    )

    # determine statepoint values and field names
    validate_file_path(args.parameters_swept, check_missing=True, check_has_extension=True, valid_extensions=('.json',))
    with args.parameters_swept.open('r') as file_params_swept:
        # params_swept = standardize_params_swept(json.load(file_params_swept))
        params_swept = json.load(file_params_swept) # NOTE: don't want to standardize here, since MixtureSpec is not a JSONSerializable type for statepoints
    
    ## identify names of statepoint and metadata fields from provided monomer dataset
    fields_not_in_df : set[str] = set(args.parameters_field).difference(monomer_df.columns)
    if any(fields_not_in_df):
        raise KeyError(f'The following requested statepoint fields not found in the provided monomer data file:\n{fields_not_in_df}')
    else: # if no extraneous fields are found, partition fields into statepoint and metadata (i.e. "other") 
        statepoint_fields = args.parameters_field
        metadata_fields = monomer_df.columns.difference(args.parameters_field)

    # create signac project directory
    project_path = (args.output_dir / args.project_name).resolve()
    validate_file_path(project_path, check_already_exists=True)
    project = init_project(args.project_name)

    ## populate data into statepoints
    for _, row in monomer_df.iterrows(): # N.B.: iterating over rows (rather than injecting into Cartesian product) since we DON'T want product along fields bundled within datafile records
        # generate job statepoints and metadata, inject shared state parameters as needed
        for param_choices in cartesian_grid(params_swept):
            job = project.open_job(
                statepoint={
                    **row.loc[statepoint_fields].to_dict(),
                    **param_choices,
                }
            )
            job.document = metadata = row.loc[metadata_fields].to_dict() # shunt accessory data from training set to document


def main() -> None:
    logging.basicConfig(level=logging.INFO, force=True)

    parser = ArgumentParser()
    parser.add_argument(
        '-mdat',
        '--monomer-data',
        type=Path,
        required=True,
        help='Path to the (formatted) monomer database file to draw records from',
    )
    parser.add_argument(
        '-num',
        '--number-to-sample',
        type=int,
        help='Optional number of records to subsample from the dataset, if the full dataset is not desired',
    )
    parser.add_argument(
        '-rand',
        '--random',
        action='store_true',
        help='When --number-to-sample is set, dictates whether the subsample should be the first N (default) or a random sample of N',
    )
    parser.add_argument(
        '-pswept',
        '--parameters-swept',
        type=Path,
        default=PARAMS_SWEPT_PATH,
        help='Path to a JSON file containing varying system parameters for project jobs',
    )
    parser.add_argument(
        '-pfield',
        '--parameters-field',
        type=str,
        nargs='+',
        default=['smiles_explicit'],
        help='List of field names from the provided monomer dataset to include as statepoint keys',
    )
    parser.add_argument(
        '-proj',
        '--project-name',
        required=True,
        help='The toplevel name of the Signac project to be initialized',
    )
    parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=_parent_dir,
        help='The directory inside which the Signac project should be initialized'
    )

    args = parser.parse_args()
    generate_statepoints(args)

if __name__ == '__main__':
    main()