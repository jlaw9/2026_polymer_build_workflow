'''For collating and storing bond length distribution data from a polymer build project'''

from rich.progress import track
from rich.logging import RichHandler

import logging
logging.basicConfig(level=logging.INFO, handlers=[RichHandler()])
LOGGER = logging.getLogger(__name__)
logging.getLogger('src.reactions').setLevel(logging.CRITICAL)

import warnings
warnings.filterwarnings(action='ignore')

from argparse import ArgumentParser, Namespace

from typing import Optional
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import Get3DDistanceMatrix

from .project import PolymerBuildProject


def extract_bond_lengths(
        project_path : Path,
        output_path : Path,
        take_first : Optional[int]=None,
    ) -> None:
    LOGGER.info(f'Loading Signac project at {project_path}')
    project = PolymerBuildProject(project_path)

    for i, job in track(enumerate(project), description='Extracting bond lengths', total=len(project)):
        if (take_first is not None) and (i > take_first):
            LOGGER.warning(f'Reached prescribed maximum of {take_first} jobs to analyze, stopping')
            break

        LOGGER.info(f'Loading molecule for job "{job.id}')
        try:
            with SDMolSupplier(job.fn(PolymerBuildProject.OLIGOMER_SDF), sanitize=False) as suppl:
                oligomer = suppl[0]
        except OSError: # catches both nonexistent and empty files - TODO: add exceptions for nonempty but malformatted molecule file?
            LOGGER.error(f'Could not read molecule from job molecule file, skipping')
            continue 
        else:
            bond_idx_pairs = tuple((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in oligomer.GetBonds())
            bond_idx_pairs = tuple(zip(*bond_idx_pairs)) # "transpose" from pairwise indices to 2 tuples or row and column indices to adhere with numpy's indexing convention
            dist_matr = Get3DDistanceMatrix(oligomer)
            bond_dists = dist_matr[bond_idx_pairs]

        with h5py.File(output_path, 'a') as hdf5_file: # need to append, otherwise only last entry is kept
            LOGGER.info(f'Caching bond length distribution')
            hdf5_file.require_dataset(job.id, shape=bond_dists.shape, dtype='f')[:] = bond_dists # allows for overwrite if field already exists
    else:
        LOGGER.info(f'Successfully extracted bond length data from project "{project_path.stem}"')

if __name__ == '__main__':
    input_parser = ArgumentParser()
    input_parser.add_argument(
        '-proj',
        '--project-path',
        type=Path,
        default=Path.cwd(),
        help='Path to the directory in which the (presumed initialized) Signac project statepoints reside',
    ),
    input_parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=None,
        help='Directory into which bond length data file should be saved (will default to project directory)'
    )
    input_parser.add_argument(
        '-of',
        '--output_file',
        type=Path,
        default='bond_lengths.hdf5',
        help='Name of the file which bond length data will be written to (must be an HDF5 file!)'
    )
    input_parser.add_argument(
        '-take',
        '--take_first',
        type=int,
        default=None,
        help='If provided as int, will take only determine bond length data for the first *this many* eligible jobs'
    )

    # post-process input arguments
    args = input_parser.parse_args()

    output_dir = args.project_path if (args.output_dir is None) else args.output_dir
    output_path = output_dir / args.output_file
    assert output_path.suffix == '.hdf5'

    extract_bond_lengths(
        project_path=args.project_path,
        output_path=output_path,
        take_first=args.take_first,
    )