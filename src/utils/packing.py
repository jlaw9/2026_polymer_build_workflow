'''Tools for generating and filling lattices of molecules'''

import numpy as np
from polymerist.maths.lattices.integral import CubicIntegerLattice


def generate_uniform_subpopulated_lattice(max_num_atoms : int, num_atoms_in_mol : int, dimension : int=3) -> CubicIntegerLattice:
    '''Create an integer lattice which accomodates a number of sites while minimizing the size of consecutive voids between empty sites'''
    num_mols = max_num_atoms // num_atoms_in_mol # NOTE: key that this is floor division and not ordinary division
    sidelen = np.ceil(num_mols**(1/dimension)).astype(int) # needed to bypass float-typing for integer-valued quantity
    sidelens = np.array([sidelen]*dimension)
    full_lattice = CubicIntegerLattice(sidelens)

    # determine how many odd and even sublattice sites to sample
    num_even_sites = full_lattice.even_idxs.size
    num_even_to_take = min(num_mols, num_even_sites)     # lower bound on occupancy in d-dims is 0.5**(d-1) (=0.25 when d=3), meaning half lattice is not guaranteed to be occupied
    num_odd_to_take  = max(0, num_mols - num_even_sites) # only choose odd sites if there are any remaining once filling the even sites

    # randomly subsample appropriate amounts of each sublattice
    even_idxs_to_keep = np.random.permutation(full_lattice.even_idxs)[:num_even_to_take] # if the even lattice is unfilled, this improves spread, and if it is full this doesn't matter
    odd_idxs_to_keep  = np.random.permutation(full_lattice.odd_idxs )[:num_odd_to_take ] # populate interstices randmoly to avoid bias towards any part of the box
    idxs_to_keep = np.concatenate([even_idxs_to_keep, odd_idxs_to_keep])
    full_lattice.points = full_lattice.points[idxs_to_keep]

    return full_lattice # TOSELF: naming here no longer makes sense as lattice is not technically full anymore: worth fixing?
