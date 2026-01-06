'''For reading, writing, and representing polymer build job statespace parameters'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
from argparse import ArgumentParser, Namespace

from typing import Any, Hashable, Mapping, TypeVar, Union
from dataclasses import dataclass, field, fields

import json
from pathlib import Path

from polymerist.genutils.fileutils.pathutils import assemble_path, allow_string_paths
from polymerist.genutils.fileutils.jsonio.jsonify import make_jsonifiable
from polymerist.genutils.fileutils.jsonio.serialize import JSONSerializable
# from polymerist.mdtools.openfftools.partialcharge.molchargers import MolCharger

from .utils.dataIO import validate_file_path
from .utils.mixtures import MixtureSpec, MixtureSpecSerializer, ParseMixtureSpec
from . import _parent_dir


# HARD-CODED PATHS WHERE PARAMETERS SHOULD LIVE
PARAMS_DIR = _parent_dir / 'parameters'
PARAMS_DIR.mkdir(exist_ok=True)
PARAMS_SWEPT_PATH  = PARAMS_DIR / 'parameters_swept.json'

# DATACLASSES TO GIVE STRUCTURE TO PARAMETER SETS
@make_jsonifiable(type_serializer=MixtureSpecSerializer)
@dataclass
class SystemParameters:
    '''
    For encapsulating a set of non-chemical parameters about a polymer build job
    e.g. related to system extent, force field configuration, mixture molecules, etc.
    '''
    DOP            : int
    n_atoms_max    : int
    pcharge_method : str = 'Espaloma-AM1-BCC' # 'NAGL'
    mixture_spec   : MixtureSpec = field(default_factory=MixtureSpec)
    # NOTE: parameters below generally shouldn't be swept through, and sensible defaults are provided for all
    forcefield              : str   = 'openff_unconstrained-2.0.0.offxml' # 'openff-2.0.0.offxml'
    minimize_oligomer       : bool  = True
    use_switching_function  : bool  = False
    switch_width_nm         : float = 0.1
    nonbonded_cutoff_nm     : float = 0.9
    box_padding_nm          : float = 0.0

def standardize_params_swept(json_dict : dict[str, JSONSerializable]) -> Mapping[str, list[Hashable]]:
    '''Read and format a JSON-serialized parameter statespace into a
    mapping from SystemParameter fields to sets of swept parameter values'''
    params_swept : dict[str, set[Any]] = dict()
    for field_ in fields(SystemParameters): # underscore to avoid confusion with dataclasses.field
        values = json_dict[field_.name] # don't use dict.get(); want a big, loud KeyError if no values for that field are provided
        if not isinstance(values, list): # reasonable assumption IFF json_dict comes directly from a JSON file
            values = [values]

        params_swept[field_.name] = [ # additional irrelevant parameters read from JSON are not injected
            field_.type(value) if not isinstance(value, field_.type) else value # assumes hashability
                for value in values
        ]
    
    return params_swept

@allow_string_paths
def write_params_json(path_config : Path, **kwargs) -> None:
    LOGGER.info('Writing swept parameters to file...')
    with path_config.open('w') as file_config:
        json.dump(
            standardize_params_swept(kwargs),
            file_config,
            indent=4,
            default=MixtureSpecSerializer.encode,
        )
    LOGGER.info(f'Swept parameters written to {path_config!s}')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, force=True)

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser')

    # parameter file path
    parser_write = subparsers.add_parser('write')
    parser_write.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=PARAMS_DIR,
        help='The directory into which the outputted data file and plots will be written'
    )
    parser_write.add_argument(
        '-namdat',
        '--name-datafile',
        type=str,
        default='parameters_swept',
    )
    parser_write.add_argument(
        '-aow',
        '--allow-overwrites',
        action='store_true',
        help='Whether to permit overwriting output files which already exist (default is False)'
    )
    # swept parameter values
    parser_write.add_argument(
        '-dop',
        '--DOP', # no long-form name to match SystemParameters signature
        type=int,
        nargs='+',
        default=[3, 5],
        help='Degree(s) of polymerization of oligomer chemistries',
    )
    parser_write.add_argument(
        '-namax',
        '--n-atoms-max',
        type=int,
        nargs='+',
        default=[20_000],
        help='Maximum number of POLYMER atoms to pack into a box',
    )
    parser_write.add_argument(
        '-pcm',
        '--pcharge-method',
        type=str,
        nargs='+',
        # choices=list(MolCharger.subclass_registry.keys()) # DEV: omitted because OpenFF import slows script to a crawl
        choices=[
            'AM1-BCC-ELF10',
            'Espaloma-AM1-BCC',
            'NAGL',
        ],
        default='Espaloma-AM1-BCC',
        help='Method to use for assigning atomic partial charges to molecules',
    )
    parser_write.add_argument(
        '-ff',
        '--forcefield',
        type=str,
        nargs='+',
        default='openff_unconstrained-2.0.0.offxml',
        # default='openff-2.0.0.offxml',
        help='Name of OpenFF ForceField to use when parameterizing molecules',
    )
    parser.add_argument(
        '-mix',
        '--mixture-spec',
        type=MixtureSpec,
        nargs='+',
        default={},
        action=ParseMixtureSpec, # bespoke processing of mixture dict
    )
    parser_write.add_argument(
        '-nemin',
        '--dont-minimize-oligomer',
        dest='minimize_oligomer', # needed to alias to 
        action='store_false', # True by default
        help='Whether to perform UFF energy minimization after generating prototype oligomer conformer'
    )
    parser_write.add_argument(
        '-usf',
        '--use_switching_function',
        action='store_true',
        help='Whether to use a switching function to smooth the non-bonded cutoff'  \
            '(will give different results between MD engines if enabled)'
    )
    parser_write.add_argument(
        '-sw-nm',
        '--switch-width-nm',
        type=float,
        nargs='+',
        default=0.1,
        help='Width of switching function interpolation (in nanometers), if switching function is enabled',
    )
    parser_write.add_argument(
        '-nbcut-nm',
        '--nonbonded-cutoff-nm',
        type=float,
        nargs='+',
        default=0.9,
        help='Distance (in nanometers) beyond which to truncate long-ranged nonbonded interactions',
    )
    parser_write.add_argument(
        '-bpad-nm',
        '--box-padding-nm',
        type=float,
        nargs='+',
        default=0.0,
        help='Extra amount (in nanometers) to pad packed box away from each face beyond tight bounding box dimensions',
    )

    # parse args from subparser and dispatch
    args = parser.parse_args()

    if args.subparser == 'write':
        ## assemble output path
        args.output_dir.mkdir(parents=False, exist_ok=True)
        path_config = assemble_path(
            args.output_dir,
            args.name_datafile,
            extension='.json',
        )
        validate_file_path(
            path_config,
            check_missing=False,
            check_already_exists=not args.allow_overwrites,
            valid_extensions=('.json',)
        )
        write_params_json(path_config, **vars(args))