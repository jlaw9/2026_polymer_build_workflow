'''For stripping out SMILES strings from assembled oligomers structures to use downstream for transfer learning'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from argparse import ArgumentParser
from rich.progress import track

from pathlib import Path
import pandas as pd

from rdkit import Chem
from polymerist.genutils.fileutils.pathutils import assemble_path
from .project import PolymerBuildProject # my custom tooling for this project

from .utils.dataIO import validate_file_path
from .project import PolymerBuildProject # my custom tooling for this project


# compile SMILES from structure files
def compile_oligomer_SMILES_data(project : PolymerBuildProject) -> pd.DataFrame:
    '''Compile SMILES strings from oligomer structure files in a PolymerBuildProject into a CSV file for ML training'''
    records = []
    for job in track(project, description='Compiling SMILES for ML training'):
        try:
            with Chem.SDMolSupplier(job.fn(PolymerBuildProject.OLIGOMER_SDF), sanitize=False, removeHs=False) as suppl:
                oligomer = suppl[0]
                Chem.SanitizeMol(oligomer)
                oligomer_H_free = Chem.RemoveHs(oligomer)
        except OSError:
            continue

        records.append({
            'hash' : job.id,
            'DOP'  : job.sp.DOP,
            'SMILES_monomer_input' : job.sp.smiles_explicit,
            'SMILES_oligomer' : Chem.MolToSmiles(oligomer, isomericSmiles=True, kekuleSmiles=False, canonical=True, allHsExplicit=False),
            'SMILES_oligomer_no_Hs' : Chem.MolToSmiles(oligomer_H_free, isomericSmiles=True, kekuleSmiles=False, canonical=True, allHsExplicit=False),
            'SMILES_oligomer_all_Hs' : Chem.MolToSmiles(oligomer, isomericSmiles=True, kekuleSmiles=False, canonical=True, allHsExplicit=True),
            'SMILES_oligomer_no_stereo' : Chem.MolToSmiles(oligomer_H_free, isomericSmiles=False, kekuleSmiles=False, canonical=True, allHsExplicit=False),
        })

    return pd.DataFrame.from_records(records)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '-path',
        '--project-path',
        type=Path,
        required=True,
        help='Path to the directory in which the (presumed initialized) Signac project statepoints reside',
    )
    parser.add_argument(
        '-od',
        '--output-dir',
        type=Path,
        default=Path.cwd(),
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
        default='oligomer_SMILES_for_ML',
    )
    args = parser.parse_args()

    # Prepare output file for write
    args.output_dir.mkdir(parents=False, exist_ok=True)
    path_smiles = assemble_path(args.output_dir, args.name_datafile, extension='.csv')
    validate_file_path(
        path_smiles,
        check_missing=False,
        check_already_exists=not args.allow_overwrites,
        valid_extensions=('.csv',),
    )

    # Compile and write SMILES data
    project = PolymerBuildProject.get_project(args.project_path)
    smiles_df = compile_oligomer_SMILES_data(project)
    
    smiles_df.to_csv(path_smiles, index=False)
