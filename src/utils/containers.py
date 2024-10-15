'''For unpacking builtin data containers'''

from typing import Any, Generator, Iterable
from itertools import product as cartesian_product


def stringify_dict(anydict : dict[Any, Any], sep : str=', ', joiner : str='-->') -> str:
    '''Create an inline string describing the mappping in a dict'''
    return sep.join(
        f'({key!s} {joiner} {value!s})'            
            for key, value in anydict.items()
    )

def cartesian_grid(param_options : dict[str, Iterable[Any]]) -> Generator[dict[str, Any], None, None]:
    '''
    Takes a dict keyed by parameter names whose values contain
    possible values for each respective parameter
    
    Exhaustively generates dicts (keyed by the same parameter names) containing every
    unique combination of those parameter values, with exactly one value for each key
    '''
    for param_point in cartesian_product(*param_options.values()):
        yield {
            param_name : param_value
                for param_name, param_value in zip(param_options.keys(), param_point)
        }