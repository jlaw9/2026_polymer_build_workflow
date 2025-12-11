from pathlib import Path
from rich.progress import track
import pandas as pd

from rdkit import Chem
from src.project import PolymerBuildProject # my custom tooling for this project


# load project from workspace
output_path = Path('oligomer_SMILES_for_ML_updated.csv')

project_dir = Path('polyID_production_LAMMPS')
# project_dir = Path('polyID_production_expanded')
project = PolymerBuildProject.get_project(project_dir)
# print(len(project))

# compile SMILES from structure files
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

df = pd.DataFrame.from_records(records)
df.to_csv(output_path, index=False)
