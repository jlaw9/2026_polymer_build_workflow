'''Initialize statepoint directories for PolyID Signac workflow'''

from argparse import ArgumentParser, Namespace

import pandas as pd
from pathlib import Path
_parent_dir = Path(__file__).parent.resolve()

from init_params import ParametersSwept, PARAMS_SWEPT_PATH 
from init_params import ParametersConfig, PARAMS_CONFIG_PATH


def read_monomer_data_fields(mono_df_path : Path) -> pd.DataFrame:
    '''Read and subselect monomer data fields according to provided read parameters'''
    ...

def generate_statepoints(
    monomer_df : pd.DataFrame,
    params_swept : ParametersSwept, 
    params_config : ParametersConfig
) -> None:
    '''Initialize project directory and job statepoint files for chosen monomer data fields'''
    ...

def init_project(args : Namespace) -> None:
    ...


def main() -> None:
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
        '-pconfig',
        '--parameters_config',
        type=Path,
        default=PARAMS_CONFIG_PATH,
        help='Path to a JSON file containing shared configuration parameters for project jobs',
    )
    parser.add_argument(
        '-pswept',
        '--parameters_swept',
        type=Path,
        default=PARAMS_SWEPT_PATH,
        help='Path to a JSON file containing varying design parameters for project jobs',
    )
    parser.add_argument(
        '-proj',
        '--project-name',
        required=True,
        help='The toplevel name of the Signac project to be initialized',
    )

    args = parser.parse_args()
    print(args)

if __name__ == '__main__':
    main()