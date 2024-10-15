'''For calculating numbers of admolecules in melts by fractions in any basis'''

from logging import getLogger
LOGGER = getLogger(__name__)

from typing import Callable, ClassVar, NewType, TypeAlias, Union
from enum import StrEnum
from dataclasses import dataclass, field

from rdkit.Chem import Descriptors, Mol, MolFromSmiles

from polymerist.smileslib.primitives import is_valid_SMILES

Smiles = NewType('Smiles', str)
Smarts = NewType('Smarts', str)
MolFract  = NewType('MolFract' , float)
MassFract = NewType('MassFract', float)
VolFract  = NewType('VolFract' , float)

AdmoleculeDict : TypeAlias = dict[Smiles, Union[MolFract, MassFract, VolFract]]

class FractionBasis(StrEnum):
    '''For choosing the calculation basis for species fraction'''
    MASS = 'mass'
    MOLAR = 'molar'
    VOLUME = 'volume' # NOTE: excluding for now, as estimating species density is pretty tricky

@dataclass
class AdmoleculeCalculator: # TODO: incorporate units?
    '''Class for stremlining and encapsulating calculations of numbers of admolecules based on fractional specification'''
    admolecules : AdmoleculeDict = field(default_factory=dict)
    basis : FractionBasis = field(default=FractionBasis.MOLAR)

    _POLYMER_KW : ClassVar[str] = '<polymer>'

    # validation methods
    def _post_init__(self) -> None:
        '''Validate input parameters'''
        assert self.admolecules_valid(self.admolecules)

    @staticmethod
    def admolecules_valid(admols : AdmoleculeDict) -> bool:
        '''Check whether an admolecule specification dict is valid'''
        admol_fract : float = 0.0 # cumulative fraction of admolecules to be added - recalculated here for portability
        for smiles, mol_fract in admols.items():
            if not is_valid_SMILES(smiles):
                LOGGER.error(f'Invalid admolecule SMILES specifier: "{smiles}"')
                return False
            if not (0 <= mol_fract <= 1):
                LOGGER.info(f'All admolecule mol fractions must be values from 0 to unity (provided {mol_fract})')
                return False
            admol_fract += mol_fract

        if not (0 <= admol_fract < 1): # must be strictly less that 1, as we also want to include some polymer component
            LOGGER.error(f'Cumulative admolecule fraction must be less than 1.0 to accomodate polymer; supplied {admol_fract}')
            return False
        
        return True
    
    @property
    def admolecule_fraction(self) -> float:
        '''The cumulative fraction of the mix taken up by all admolecules'''
        return sum(self.admolecules.values())

    ## aliases
    @property
    def admols(self) -> AdmoleculeDict:
        '''Alias of self.admolecules for convenience'''
        return self.admolecules

    @property
    def admol_fract(self) -> float:
        '''Alias of self.admolecule_fraction for convenience'''
        return self.admolecule_fraction
    
    # calculation methods
    def _num_molecules_molar(self, num_polymers : int) -> dict[Smiles, int]:
        '''Calculate the number of admolecules and polymers needed on a molar basis'''
        polymer_mol_fract : float = 1.0 - self.admol_fract
        N_total : int = round(num_polymers / polymer_mol_fract)

        num_mols = {
            smiles : round(N_total * mol_fract)
                for smiles, mol_fract in self.admols.items()
        }
        num_mols[self._POLYMER_KW] = num_polymers

        return num_mols

    def _num_molecules_mass(self, num_polymers : int, polymer_mass : float) -> dict[Smiles, int]:
        '''Calculate the number of admolecules and polymers needed on a molar basis'''
        polymer_mass_fract : float = 1.0 - self.admol_fract
        polymer_mass : float = num_polymers * polymer_mass
        M_total : float = polymer_mass / polymer_mass_fract
        
        num_mols : dict[Smiles, float] = {}
        for smiles, mass_fract in self.admols.items():
            Mw = Descriptors.ExactMolWt(MolFromSmiles(smiles)) # moleular weight of the current species
            num_mols[smiles] = round( (M_total * mass_fract)/Mw )
        num_mols[self._POLYMER_KW] = num_polymers

        return num_mols
    
    def _num_molecules_volume(self, num_polymers : int, polymer_density : float) -> dict[Smiles, int]:
        '''Calculate the number of admolecules and polymers needed on a molar basis'''
        raise NotImplementedError
    
    _CALC_METHODS : ClassVar[dict[FractionBasis, Callable[['AdmoleculeCalculator', int, Smiles], dict[Smiles, int]]]] = {
        FractionBasis.MASS  : _num_molecules_mass,
        FractionBasis.MOLAR : _num_molecules_molar,
        FractionBasis.VOLUME : _num_molecules_volume,
    }

    def num_molecules(self, num_polymers : int, *args, **kwargs) -> dict[Smiles, int]:
        '''
        Calculate the number of each admolecule needed to meet the given admolecule fractions
        as closely as possible for the chosen basis given the target number of (assumed-to-be
        identical) polymers in the melt and the polymer SMILES
        '''
        calc_method = self._CALC_METHODS[self.basis]
        return calc_method(self, num_polymers, *args, **kwargs)