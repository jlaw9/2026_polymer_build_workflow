'''For compiling data about ring piercing from a polymer structure build project'''

from rich.progress import track
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import Get3DDistanceMatrix

from src.project import PolymerBuildProject
from src.utils.piercing import assess_ring_piercing


# plotting params
nrows : int = 1
ncols : int = 2
aspect : float = 5/4
scale : float = 6.0

fontsize : int = 14
n_bins : int = 75

ANGSTROM : str = "\u212B"

# load project
project_path = Path('polyID_production')
project = PolymerBuildProject.get_project(project_path)
print(len(project))

# compile piercing data
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

    ring_piercing_idxs : dict[tuple[int], tuple[tuple[int, int]]] = assess_ring_piercing(oligomer)
    n_rings_pierced = sum(bool(idxs) for idxs in ring_piercing_idxs.values())
    found_pierced_rings = any(ring_piercing_idxs.values())

    records.append({
        'job_hash' : job.id,
        'bond_dist_min' : bond_dists.min(),
        'bond_dist_max' : bond_dists.max(),
        'found_pierced_rings' : found_pierced_rings,
        'n_rings_pierced' : n_rings_pierced,
    })

ring_piercing_df = pd.DataFrame.from_records(records)
pierced_subset = ring_piercing_df[ring_piercing_df['found_pierced_rings']]
pierced_subset.to_csv('pierced_systems.csv')
print(f'{len(ring_piercing_df)} ring-pierced systems found')

# produce plots, save data
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

fig.savefig('bond_lengths_pierced_extremal.png', bbox_inches='tight')
plt.close()