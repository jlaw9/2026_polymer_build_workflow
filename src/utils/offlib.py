'''For reading info from and writing info to OpenFF Molecule/Topology objects'''

from typing import Union

import re
from collections import Counter

from openff.toolkit import Molecule, Topology
from openff.units.elements import SYMBOLS


HILL_REGEX = re.compile('(?P<element>[A-Z][a-z]?)(?P<count>[0-9]*)') # break apart Hill formula into just unique elements (one capital letter, one or no lowercase letters, any (including none) digits)
def elem_counts_hill(offmol : Molecule) -> dict[str, int]:
    '''Extract unique elements and their counts from a Molecule object's Hill formula'''
    elem_counts = {}
    for match in re.finditer(HILL_REGEX, offmol.to_hill_formula()):
        match_groups = match.groupdict()
        if match_groups['count'] == '':
            match_groups['count'] = '1'
        elem_counts[match_groups['element']] = int(match_groups['count'])
    return elem_counts

def elem_counts(offobj : Union[Molecule, Topology], from_hill_formula : bool=False) -> dict[str, int]:
    '''Takes an penFF Molecule or Topology object and returns a dict keyed by
    unique element symbols whose values count the number of occurrences of that element'''

    if from_hill_formula:
        if isinstance(offobj, Topology):
            raise ValueError(f'Cannot parse Hill formula from {Topology!s} object which does not implement to_hill_formula()')
        # logging.warn('Extracting atom counts from Hill formula, rather than direct count')
        return elem_counts_hill(offobj)
    else:
        elem_counts = Counter(
            SYMBOLS[atom.atomic_number]
                for atom in offobj.atoms # this works because both Molecule and Topology implement the "atoms" iterator
        )
        return dict(elem_counts)