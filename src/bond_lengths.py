'''For collating and storing bond length distribution data from a polymer build project'''

from rich.progress import track
from rich.logging import RichHandler

import logging
logging.basicConfig(level=logging.INFO, handlers=[RichHandler()])
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

import warnings
warnings.filterwarnings(action='ignore')

from argparse import ArgumentParser, Namespace

from typing import Sequence, Optional
from pathlib import Path

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import Get3DDistanceMatrix

from project import PolymerBuildProject
from polymerist.genutils.fileutils.pathutils import assemble_path

logging.getLogger('numexpr.utils').setLevel(logging.CRITICAL) # suppress numpy logs
logging.getLogger('reactions').setLevel(logging.CRITICAL) # suppress logs from src.reactions (called on .project import)


DEFAULT_DATAFILE_NAME : str = 'bond_lengths'

ANGSTROM : str = "\u212B"
LEQ : str = '<='

## taken from https://sites.google.com/site/chempendix/bond-lengths
EXPER_BOND_LENGTHS : dict[str, float] = { 
    'O-H' : 0.97,
    'C-H' : 1.09,
    'C=O' : 1.23,
    'C=C' : 1.34,
    # 'C:C' : 1.39,
    'C-N=' : 1.4,
    'C-O' : 1.43,
    'C-C' : 1.54,
    'C=S' : 1.73,
    'C-S' : 1.83,
}

def extract_bond_lengths(
    project_path : Path,
    output_path : Path,
    take_first_n : Optional[int]=None,
) -> None:
    '''Compile data on all lengths of all bonds of all polymers in a given project and store in HDF5 file'''
    LOGGER.info(f'Loading Signac project at {project_path}')
    project = PolymerBuildProject(project_path)

    for i, job in track(
        enumerate(project),
        description='Extracting bond lengths',
        total=len(project) if (take_first_n is None) else take_first_n,
    ):
        if (take_first_n is not None) and (i > take_first_n):
            LOGGER.warning(f'Reached prescribed maximum of {take_first_n} jobs to analyze, stopping')
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
            bond_idx_pairs = tuple(zip(*bond_idx_pairs)) # "transpose" from pairwise indices to 2 tuples of row and column indices to adhere with numpy's indexing convention
            dist_matr = Get3DDistanceMatrix(oligomer)
            bond_dists = dist_matr[bond_idx_pairs]

        with h5py.File(output_path, 'a') as hdf5_file: # need to append, otherwise only last entry is kept
            LOGGER.info(f'Caching bond length distribution')
            hdf5_file.require_dataset(job.id, shape=bond_dists.shape, dtype='f')[:] = bond_dists # allows for overwrite if field already exists
    else:
        LOGGER.info(f'Successfully extracted bond length data from project "{project_path.stem}"')

def plot_bond_length_distribution(
    all_bond_dists : np.ndarray[float],
    ticks_per_angstrom_binned : Optional[Sequence[int]]=None,
    bond_bin_edges : Optional[Sequence[float]]=None,
    number_of_rows : int=1,
    scale : float=7.0,
    aspect : float=1.2,
    fontsize : float = 15.0,
    plot_experimental : bool=False,
) -> tuple[Figure, Axes]:
    '''Plot bond length distribution from cached HDF5 data file'''
    # Determine partition of bond length range to cut histogram into
    if ticks_per_angstrom_binned is None:
        ticks_per_angstrom_binned = [10, 25, 5]

    if bond_bin_edges is None:
        bond_bin_edges = [0.9, 1.6]
    bond_bin_edges = np.array(bond_bin_edges)

    assert len(ticks_per_angstrom_binned) == (len(bond_bin_edges) + 1)
    bond_bin_idxs = np.digitize(all_bond_dists, bins=bond_bin_edges)
    bond_bin_labels = np.unique(bond_bin_idxs)

    number_of_columns = np.max([len(bond_bin_labels), len(bond_bin_edges) + 1])
    fig, axes = plt.subplots(number_of_rows, number_of_columns, figsize=(number_of_columns*scale*aspect, number_of_rows*scale))

    # Histogram bond lengths
    for i, bin_label in enumerate(bond_bin_labels):
        ax = axes[i]
        ticks_per_angstrom = ticks_per_angstrom_binned[i]
        
        # determine number of histogram bins (distinct from range bins!) by square root rule
        bond_dists_in_bin = all_bond_dists[bond_bin_idxs == bin_label]
        n_hist_bins : int = np.floor(np.sqrt(len(bond_dists_in_bin))).astype(int)
        
        lower = '' if (i == 0) else f'{bond_bin_edges[i - 1]} {LEQ} '
        upper = '' if (i >= len(bond_bin_edges)) else f' {LEQ} {bond_bin_edges[i]}'
        bin_desc = f'{lower}d{upper}'

        # set x-axis ticks based on range and specified number
        tick_idx_min = np.floor(bond_dists_in_bin.min()*ticks_per_angstrom)
        tick_idx_max = np.ceil(bond_dists_in_bin.max()*ticks_per_angstrom)
        xticks = np.linspace(
            tick_idx_min/ticks_per_angstrom,
            tick_idx_max/ticks_per_angstrom,
            num=int(tick_idx_max - tick_idx_min + 1),
        )
        # construct histogram for bin
        ax.hist(bond_dists_in_bin, bins=n_hist_bins)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f'{tick:1.2f}' for tick in xticks], rotation=-30)

        ax.set_title(f'Bond length distribution ({bin_desc})', fontsize=fontsize)
        ax.set_xlabel(f'Bond length ({ANGSTROM})', fontsize=0.8*fontsize)
        ax.set_ylabel('Number of bonds', fontsize=fontsize)

    if plot_experimental:
        # plot experimental bond lengths
        cmap = plt.get_cmap('turbo')
        diatomic_bond_colors = iter(cmap(np.linspace(0, 1, num=len(EXPER_BOND_LENGTHS))))

        for pairtype, diatomic_bond_length in EXPER_BOND_LENGTHS.items():
            ax = axes[np.digitize(diatomic_bond_length, bins=bond_bin_edges)]
            ax.axvline(diatomic_bond_length, color=next(diatomic_bond_colors), linestyle='--', label=pairtype)
            _ = ax.legend(fontsize=fontsize)

    return fig, axes


if __name__ == '__main__':
    # CLI args
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser') # need to keep track (rather than using default arg.func) since input are also post-processed differently

    ## Bond length extraction
    extract_parser = subparsers.add_parser('extract', help='Extract and cache bond length distribution data for a polymer project')
    extract_parser.add_argument(
        '-path',
        '--project-path',
        type=Path,
        default=Path.cwd(),
        help='Path to the directory in which the (presumed initialized) Signac project statepoints reside',
    )
    extract_parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=None,
        help='Directory into which bond length data file should be saved (will default to project directory)'
    )
    extract_parser.add_argument(
        '-on',
        '--output-name',
        type=Path,
        default=DEFAULT_DATAFILE_NAME,
        help='Name of the file which bond length data will be written to (must be an HDF5 file!)'
    )
    extract_parser.add_argument(
        '-take',
        '--take-first-n',
        type=int,
        default=None,
        help='If provided as int, will take only determine bond length data for the first *this many* eligible jobs'
    )

    ## Bond length plotting 
    plot_parser = subparsers.add_parser('plot', help='Plot bond length distribution from cached data file')
    plot_parser.add_argument(
        '-dp',
        '--data-path',
        type=Path,
        default=assemble_path(Path.cwd(), f'{DEFAULT_DATAFILE_NAME}', extension='hdf5'),
    )
    plot_parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='Directory into which bond length distribution plots should be saved (will default to current working directory)'
    )
    plot_parser.add_argument(
        '-on',
        '--output-name',
        type=Path,
        default='bond_length_distribution',
        help='Name of the file which bond length distribution plots will be saved as (will be saved as PNG files)'
    )
    plot_parser.add_argument(
        '-binedg',
        '--bond-bin-edges',
        type=float,
        nargs='+',
        default=[0.9, 1.6],
        help='Bond lengths (in Angstrom) at which to split up the bond length distribution into separate subranges',
    )
    plot_parser.add_argument(
        '-tpang',
        '--ticks-per-angstrom',
        type=int,
        nargs='+',
        default=[10, 25, 5],
        help='Number of ticks per Angstrom to render per plot on each respective binned subrange of bond lengths'
    )
    plot_parser.add_argument(
        '-nrow',
        '--number-of-rows',
        type=int,
        default=1,
        help='Number of rows in which to arrange subplots',
    )
    plot_parser.add_argument(
        '-s',
        '--scale',
        type=float,
        default=7.0,
        help='Scaling factor for figure size',
    )
    plot_parser.add_argument(
        '-asp',
        '--aspect',
        type=float,
        default=1.2,
        help='Aspect ratio for figure size',
    )
    plot_parser.add_argument(
        '-fs',
        '--fontsize',
        type=float,
        default=15.0,
        help='Font size for plot text elements',
    )
    plot_parser.add_argument(
        '-exper',
        '--plot-experimental-bond-lengths',
        action='store_true',
        help='Whether to overlay experimental diatomic bond lengths for reference',
    )

    # post-process input arguments
    args = parser.parse_args()
    if args.subparser == 'extract':
        output_dir = args.project_path if (args.output_dir is None) else args.output_dir
        extract_bond_lengths(
            project_path=args.project_path,
            output_path=assemble_path(output_dir, args.output_name, extension='hdf5'),
            take_first_n=args.take_first_n,
        )
    elif args.subparser == 'plot':
        # read bond length data
        all_bond_dists_sequential : list[np.ndarray[float]] = []
        with h5py.File(args.data_path, 'r') as hdf5_file:
            for indiv_bond_dists in hdf5_file.values():
                all_bond_dists_sequential.append(indiv_bond_dists[:])
        all_bond_dists : np.ndarray[float] = np.concatenate(all_bond_dists_sequential)

        # plot and save
        plot_path = assemble_path(args.output_dir, args.output_name, extension='png')
        fig, ax = plot_bond_length_distribution(
            all_bond_dists=all_bond_dists,
            ticks_per_angstrom_binned=args.ticks_per_angstrom,
            bond_bin_edges=args.bond_bin_edges,
            number_of_rows=args.number_of_rows,
            scale=args.scale,
            aspect=args.aspect,
            fontsize=args.fontsize,
            plot_experimental=args.plot_experimental_bond_lengths,
        )
        fig.savefig(plot_path, bbox_inches='tight')
