'''Unit tests for smilesparse.py'''

from typing import Sequence
import pytest

from ..cheminf import parse_monomer_smiles


monomer_smiles_valid = [ # all of the following format variations will be accepted
    'O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1',
    ('O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'),
    "('O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1')",
    ['O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'],
    "['O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1']",
]
monomer_smiles_invalid = [ # all of the following format variations will be accepted
    'O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1',
    ('O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'),
    "('O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1')",
    ['O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'],
    "['O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1']",
]

@pytest.mark.parametrize('valid_smiles', monomer_smiles_valid)
def test_parse_smiles_valid(valid_smiles : Sequence[str]) -> None:
    '''Test that incorrectly-formatted SMILES fields are processed'''
    assert parse_monomer_smiles(valid_smiles) is not None

@pytest.mark.parametrize('invalid_smiles', monomer_smiles_invalid)
def test_parse_smiles_invalid(invalid_smiles : Sequence[str]) -> None:
    '''Test that incorrectly-formatted SMILES fields are rejected'''
    assert parse_monomer_smiles(invalid_smiles) is None