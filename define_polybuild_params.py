'''Programmatically redefine swept and shared parameter sets for a polymer build Signac Project'''

import json
from pathlib import Path


PARAMS_SWEPT_PATH = Path('parameters_swept.json')
PARAMS_SHARED_PATH = Path('parameters_shared.json')


params_swept = { # define other parameters to sweep here!
    'DOP' : [
        3,
        # 5,
    ],
    'n_atoms_max' : [
        10_000,
        # 20_000,
    ],
    'pcharge_method' : [
        'Espaloma-AM1-BCC',
        # 'NAGL',
    ],
}
with PARAMS_SWEPT_PATH.open('w') as file:
    json.dump(params_swept, file, indent=4)

params_shared = { # shared default parameters that we don't expect to have to sweep through
    'minimize_oligomer'      : True,
    'forcefield'            : 'openff_unconstrained-2.0.0.offxml', # 'openff-2.0.0.offxml',
    'use_switching_function' : False,
    'switch_width_nm'        : 0.1,
    'nonbonded_cutoff_nm'    : 0.9,
    'box_padding_nm'         : 0.0,
}
with PARAMS_SHARED_PATH.open('w') as file:
    json.dump(params_shared, file, indent=4)