'''Unit tests for smilesparse.py'''

from format_monomer_data import parse_monomer_smiles

def test_parse_smiles_invalid() -> None:
    '''Test that incorrectly-formatted SMILES fields are processed'''
    mono_smiles_accepted = [ # all of the following format variations will be accepted
        'O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1',
        ('O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'),
        "('O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1')",
        ['O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'],
        "['O=C(Cl)Cl.Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1']",
    ]

    for smi in mono_smiles_accepted:
        assert(parse_monomer_smiles(smi) is not None)

def test_parse_smiles_invalid() -> None:
    '''Test that incorrectly-formatted SMILES fields are rejected'''
    mono_smiles_invalid = [ # all of the following format variations will be accepted
        'O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1',
        ('O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'),
        "('O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1')",
        ['O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1'],
        "['O=C(Cl)Cl, Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1']",
    ]
    
    for smi in mono_smiles_invalid:
        assert(parse_monomer_smiles(smi) is None)