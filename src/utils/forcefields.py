'''Helper methods for reading and writing OpenFF SMIRNOFF force field files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable, Union

import re
from pathlib import Path
from functools import reduce

from openff.toolkit.utils.exceptions import SMIRNOFFParseError
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField, _get_installed_offxml_dir_paths

from polymerist.genutils.fileutils.pathutils import allow_pathlib_paths


FF_KEYWORD : str = '$FORCEFIELDS'
FF_DIRNAME : str = 'forcefields'


LOCAL_FF_PATTERN = re.compile(fr'{re.escape(FF_KEYWORD)}/(.*)')

@allow_pathlib_paths
def _substitute_ff_dir_path(ff_path : str, local_ff_dir : Path) -> str:
    '''If the forcefield keyword is found at the start of a path, replace it with the absolute path to the local forcefields directory'''
    ffmatch = re.match(LOCAL_FF_PATTERN, ff_path) # N.B.: opted for match to ONLY hit on start of string (keyword in middle doesn't count!)
    if ffmatch is None:
        return ff_path
    else:
        path_remainder, = ffmatch.groups()
        return str(local_ff_dir / path_remainder)
    
def load_composite_forcefield(
    *sources : Iterable[Union[str, Path]],
    local_ff_dir : Path,
    deregister_am1bcc : bool=False,    
) -> ForceField:
    '''
    Load an OpenFF Forcefield from mixed SMIRNOFF-spec file sources
    Intelligently processes locally-defined forcefields regardless of CWD
    '''
    forcefield = reduce(
        ForceField.combine,
        (ForceField(_substitute_ff_dir_path(ffpath, local_ff_dir=local_ff_dir)) for ffpath in sources),
        ForceField(), # init with enpty forcefield as base case
    )
    if deregister_am1bcc and ('ToolkitAM1BCC' in forcefield.registered_parameter_handlers):
        forcefield.deregister_parameter_handler('ToolkitAM1BCC') # forcibly remove AM1BCC handler so a failed isomorphism doesn't result in prohibitively-long AM1BCC calculation
    
    return forcefield