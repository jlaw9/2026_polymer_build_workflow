'''Programmatically redefine swept and shared parameter sets for a polymer build Signac Project'''

import logging
logging.basicConfig(level=logging.INFO)

from dataclasses import dataclass, field

from polymerist.genutils.fileutils.jsonio.jsonify import make_jsonifiable
from . import _parent_dir


# HARD-CODED PATHS WHERE PARAMETERS SHOULD LIVE
PARAMS_DIR = _parent_dir / 'parameters'
PARAMS_DIR.mkdir(exist_ok=True)

PARAMS_CONFIG_PATH = PARAMS_DIR / 'parameters_config.json'
PARAMS_SWEPT_PATH  = PARAMS_DIR / 'parameters_swept.json'

# DATACLASSES TO GIVE STRUCTURE TO PARAMETER SETS
@make_jsonifiable
@dataclass
class ParametersSwept:
    '''Encapsulation class for tracking varying design parameters between polymers'''
    DOP             : list[int] = field(default_factory=list)
    n_atoms_max     : list[int] = field(default_factory=list)
    pcharge_method  : list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        for attrname in ('DOP', 'n_atoms_max', 'pcharge_method'):
            if not getattr(self, attrname):
                raise ValueError(f'Required attribute "{attrname}" unset')

@make_jsonifiable
@dataclass
class ParametersConfig:
    '''Encapsulation class for tracking fixed configuration parameters that shouldn't be changed'''
    forcefield              : str
    minimize_oligomer       : bool
    use_switching_function  : bool = False
    switch_width_nm         : float = 0.1
    nonbonded_cutoff_nm     : float = 0.9
    box_padding_nm          : float = 0.0


# PROCEDURAL REGENERATION OF PARAMETER FILES IF THIS SCRIPT IS INVOKED DIRECTLY
if __name__ == '__main__':
    params_swept = ParametersSwept( # define other parameters to sweep here!
        DOP=[
            3,
            # 5,
        ],
        n_atoms_max=[
            10_000,
            # 20_000,
        ],
        pcharge_method=[
            'Espaloma-AM1-BCC',
            # 'NAGL',
        ],
    )
    logging.info('Writing swept parameters to file')
    params_swept.to_file(PARAMS_SWEPT_PATH)

    params_config = ParametersConfig( # shared default parameters that we don't expect to have to sweep through
        forcefield='openff_unconstrained-2.0.0.offxml', # 'openff-2.0.0.offxml',
        minimize_oligomer=True,
        use_switching_function=False,
        switch_width_nm=0.1,
        nonbonded_cutoff_nm=0.9,
        box_padding_nm=0.0,
    )
    logging.info('Writing config parameters to file')
    params_config.to_file(PARAMS_CONFIG_PATH)