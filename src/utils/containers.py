'''For unpacking builtin data containers'''

from typing import Any


def stringify_dict(anydict : dict[Any, Any], sep : str=', ', joiner : str='-->') -> str:
    '''Create an inline string describing the mappping in a dict'''
    return sep.join(
        f'({key!s} {joiner} {value!s})'            
            for key, value in anydict.items()
    )