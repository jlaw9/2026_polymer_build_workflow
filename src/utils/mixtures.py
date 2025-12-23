'''For representing mixtures of molecules, and converting between common bases for amount of species'''

from typing import Any, Hashable, Optional
from collections import UserDict
from polymerist.genutils.fileutils.jsonio.serialize import TypeSerializer


class MixtureSpec(UserDict):
    '''Hashable, immutable mapping from species SMILES string to amount of that species in a chosen basis'''
    def __init__(self, initdict : Optional[dict[Hashable, Any]]=None, /,  **kwargs) -> None:
        self._locked = False

        if initdict is None:
            super().__init__(**kwargs)
        else:
            initdict = {
                Chem.CanonSmiles(smi) : amount # TODO: eventually, add support for choice of basis (as (amount, basis) pairs)
                    for smi, amount in initdict.items()
            }
            super().__init__(initdict, **kwargs)

        self._locked = True

    def __setitem__(self, key : Hashable, value : tuple[float]) -> None:
        if self._locked:
            raise PermissionError(f"Direct key-value assignment is not allowed")
        else:
            super().__setitem__(key, value)

    def __hash__(self) -> str:
        return hash(tuple(sorted(
            (smi, amount)
                for smi, amount in self.items()
        )))

class MixtureSpecSerializer(TypeSerializer, python_type=MixtureSpec):
    '''For gracefully serializing MixtureSpecs against JSON-formatted data'''
    @staticmethod
    def encode(python_obj : MixtureSpec) -> dict:
        return dict(python_obj)

    @staticmethod
    def decode(json_obj : dict) -> MixtureSpec:
        return MixtureSpec(json_obj)