'''For compiling data about ring piercing from a polymer structure build project'''

from argparse import ArgumentParser
from rich.progress import track

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from rdkit import Chem
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import Get3DDistanceMatrix

from project import PolymerBuildProject
from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.rdutils.rdcoords.piercing import summarize_ring_piercing


ANGSTROM : str = "\u212B"

def collate_project_ring_piercing_data(project : PolymerBuildProject) -> pd.DataFrame:
    '''
    Takes a PolymerBuildProject and compiles a DataFrame summarizing
    bond length distribution and ring piercing detection data
    across all jobs within the project
    '''
    records = []
    for job in track(project, description='Detecting ring piercing'):
        try:
            with SDMolSupplier(job.fn(PolymerBuildProject.OLIGOMER_SDF), sanitize=False, removeHs=False) as suppl:
                oligomer = suppl[0]
                Chem.SanitizeMol(oligomer)
        except OSError:
            continue

        bond_idx_pairs = tuple((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in oligomer.GetBonds())
        bond_idx_pairs = tuple(zip(*bond_idx_pairs)) # "transpose" from pairwise indices to 2 tuples of row and column indices to adhere with numpy's indexing convention
        dist_matr = Get3DDistanceMatrix(oligomer)
        bond_dists = dist_matr[bond_idx_pairs]

        ring_piercing_idxs : dict[tuple[int], tuple[tuple[int, int]]] = summarize_ring_piercing(oligomer)
        n_rings_pierced = sum(bool(idxs) for idxs in ring_piercing_idxs.values())
        found_pierced_rings = any(ring_piercing_idxs.values())

        records.append({
            'job_hash' : job.id,
            'bond_dist_min' : bond_dists.min(),
            'bond_dist_max' : bond_dists.max(),
            'found_pierced_rings' : found_pierced_rings,
            'n_rings_pierced' : n_rings_pierced,
        })

    return pd.DataFrame.from_records(records)

def plot_pierced_bond_length_distribution(
    pierced_subset : pd.DataFrame,
    nrows : int=1,
    ncols : int=2,
    scale : float=6.0,
    aspect : float=5/4,
    fontsize : int=14,
    n_bins : int=75,
) -> tuple[Figure, Axes]:
    '''Plot distribution of bond lengths in ring-pierced systems'''
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols/nrows*aspect*scale, scale))
    fig.suptitle(f'Extremal bond distances in ring-pierced systems ({len(pierced_subset)} pierced oligomers detected)')

    plot_info = {
        'Minimal Bond Distance' : ('bond_dist_min',),
        'Maximal Bond Distance' : ('bond_dist_max',),
    }

    for ax, (subplot_title, (data_field,)) in zip(axes.flatten(), plot_info.items()):
        ax.hist(pierced_subset[data_field], bins=n_bins)

        ax.set_title(subplot_title)
        ax.set_xlabel(f'Bond length ({ANGSTROM})', fontsize=0.8*fontsize)
        ax.set_ylabel('#systems with bond length', fontsize=fontsize)

    return fig, axes


if __name__ == '__main__':
    # Solicit inputs from CLI
    parser = ArgumentParser()
    parser.add_argument(
        '-path',
        '--project-path',
        type=Path,
        required=True,
        help='Path to the directory in which the (presumed initialized) Signac project statepoints reside',
    )
    parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='The directory into which the outputted data file and plots will be written'
    )
    parser.add_argument(
        '-aow',
        '--allow-overwrites',
        action='store_true',
        help='Whether to permit overwriting output files which already exist (default is False)'
    )
    parser.add_argument(
        '-namdat',
        '--name-datafile',
        type=str,
        default='ring_piercing_data',
    )
    parser.add_argument(
        '-nplot',
        '--name-plot',
        type=str,
        default='ring_piercing_plot',
    )
    parser.add_argument(
        '-nrows',
        '--number-of-rows',
        type=int,
        default=1,
        help='Number of rows in the output plot',
    )
    parser.add_argument(
        '-ncols',
        '--number-of-columns',
        type=int,
        default=2,
        help='Number of columns in the output plot',
    )
    parser.add_argument(
        '-sc',
        '--scale',
        type=float,
        default=6.0,
        help='Scaling factor for the output plot size',
    )
    parser.add_argument(
        '-asp',
        '--aspect',
        type=float,
        default=5/4,
        help='Aspect ratio (width/height) of the output plot',
    )
    parser.add_argument(
        '-fs',
        '--fontsize',
        type=int,
        default=14,
        help='Base font size for the output plot',
    )
    parser.add_argument(
        '-nb',
        '--number-of-bins',
        type=int,
        default=75,
        help='Number of bins to use in the bond length histograms',
    )
    
    args = parser.parse_args()
    args.output_dir.mkdir(parents=False, exist_ok=args.allow_overwrites)

    # Extract data
    project = PolymerBuildProject.get_project(args.project_path)
    ring_piercing_df = collate_project_ring_piercing_data(project)
    pierced_subset = ring_piercing_df[ring_piercing_df['found_pierced_rings']]
    
    path_piercing_data = assemble_path(args.output_dir, args.name_datafile, extension='.csv')
    pierced_subset.to_csv(path_piercing_data)

    # Plot bond length distribution
    fig, ax = plot_pierced_bond_length_distribution(
        pierced_subset,
        nrows=args.number_of_rows,
        ncols=args.number_of_columns,
        scale=args.scale,
        aspect=args.aspect,
        fontsize=args.fontsize,
        n_bins=args.number_of_bins,
    )
    path_plot = assemble_path(args.output_dir, args.name_plot, extension='.png')
    fig.savefig(path_plot, bbox_inches='tight')
    plt.close()