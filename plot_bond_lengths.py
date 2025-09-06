'''Script for generating easier-to-parse visualization of bond length distributions for a project'''

import h5py
from pathlib import Path

from math import sqrt, floor
import numpy as np
import matplotlib.pyplot as plt


# specify whicch data we're looking at
bond_length_data : Path = Path('polyID_production/bond_lengths.hdf5')

# set plot parameters
scale : float = 7.0
aspect : float = 1.2
ticks_per_angstrom_binned : list[int] = [10, 25, 5]
bond_length_bins = np.array([0.9, 1.6]) # separators along which to partition domain of bond lengths for plotting
fontsize : float = 15.0

# read in bond length distribution data
all_bond_dists : list[np.ndarray[float]] = []
with h5py.File(bond_length_data, 'r') as hdf5_file:
    for indiv_bond_dists in hdf5_file.values():
        all_bond_dists.append(indiv_bond_dists[:])
all_bond_dists : np.ndarray[float] = np.concatenate(all_bond_dists)

# partition distrubtion into bins
ANGSTROM : str = "\u212B"
bond_bin_idxs = np.digitize(all_bond_dists, bins=bond_length_bins)
bond_bin_labels = np.unique(bond_bin_idxs)

# create blank canvas
n_rows : int = 1
n_cols = len(bond_bin_labels)
fig, axes = plt.subplots(1, n_cols, figsize=(n_cols*scale*aspect, n_rows*scale))
for i, bin_label in enumerate(bond_bin_labels):
    # extract lengths falling in the target range
    ax = axes[i]
    ticks_per_angstrom = ticks_per_angstrom_binned[i]
    bond_dists_in_bin = all_bond_dists[bond_bin_idxs == bin_label]
    n_hist_bins = floor(sqrt(len(bond_dists_in_bin))) # determine number of histogram bins (distinct from range bins!) by square root rule
    
    lower = '' if (i == 0) else f'{bond_length_bins[i - 1]} <= '
    upper = '' if (i >= len(bond_length_bins)) else f' <= {bond_length_bins[i]}'
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

# plot experimental bond lengths
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
cmap = plt.get_cmap('turbo')
diatomic_bond_colors = iter(cmap(np.linspace(0, 1, num=len(EXPER_BOND_LENGTHS))))

for pairtype, diatomic_bond_length in EXPER_BOND_LENGTHS.items():
    ax = axes[np.digitize(diatomic_bond_length, bins=bond_length_bins)]
    ax.axvline(diatomic_bond_length, color=next(diatomic_bond_colors), linestyle='--', label=pairtype)
    _ = ax.legend(fontsize=fontsize)

# plt.show()

fig.savefig('bond_length_distribution.png', bbox_inches='tight')