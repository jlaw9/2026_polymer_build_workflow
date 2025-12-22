from argparse import ArgumentParser
from pathlib import Path

from dataclasses import dataclass
import pandas as pd

from polymerist.genutils.fileutils.pathutils import assemble_path


@dataclass
class MonomerInfo:
    name : str
    monomer_smiles : tuple[str, ...]
    acronym : str = ''

preparations = [
    MonomerInfo(
        'polyethylene',
        ('C=C',),
        'PE',
    ),
    MonomerInfo(
        'poly(vinylidene-dichloride)',
        ('C=C(-Cl)-Cl',),
        'PVDC',
    ),
    MonomerInfo( # also ambiguous; how to get 6:1 ratio??
        'poly(ethylene vinyl alcohol)',
        ('C=C', 'CC(=O)OC=C'),
        'EVOH',
    ),
    MonomerInfo(
        'poly(ethylene terephthalate)',
        ('COC(=O)c1ccc(cc1)C(=O)OC', 'OCCO'),
        'PET',
    ),
    MonomerInfo(
        'poly(ethylene furanoate)',
        ('OC(=O)c1ccc(o1)C(=O)O', 'OCCO'),
        'PEF',
    ),
    MonomerInfo(
        'poly(L-lactic acid)',
        ('C[C@@H](C(=O)O)O',),
        'PLLA',
    ),
    MonomerInfo(
        'polystyrene',
        ('c1ccccc1-C=C',),
        'PS',
    ),
    MonomerInfo(
        'poly(acrylonitrile)',
        ('N#CC=C',),
        'PAN',
    ),
    MonomerInfo(
        'Nylon 6',
        ('[NH2]CCCCCC(=O)O',),
        'PA6',
    ),
    MonomerInfo(
        'poly(3-hydroxybutyrate)',
        ('OC(C)CC(=O)O',),
        'P3HB',
    ),
    MonomerInfo(
        'poly(glycolic acid)',
        ('OCC(=O)O',),
        'PGA',
    ),
    MonomerInfo( # ambiguous
        'poly(lactic-co-glycolic acid)',
        ('OC(C)C(=O)O', 'OCC(=O)O'),
        'PLGA',
    ),
]

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=Path('monomer_data_raw'),
        help='The directory into which the outputted data file and plots will be written'
    )
    parser.add_argument(
        '-aow',
        '--allow-overwrites',
        action='store_true',
        help='Whether to permit overwriting output files which already exist (default is False)'
    )
    parser.add_argument(
        '-namdat',
        '--name-datafile',
        type=str,
        default='permeability_pilot_set',
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=False, exist_ok=args.allow_overwrites)

    perm_pilot_df = pd.DataFrame.from_records([vars(monoinfo) for monoinfo in preparations])
    perm_pilot_path = assemble_path(args.output_dir, args.name_datafile, extension='.csv')
    perm_pilot_df.to_csv(perm_pilot_path, index=False)