import pandas as pd
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 300,300

class Molecule:
    def __init__(self, mol: Type[rdkit.Chem.Mol] = None, smiles: str = None) -> None:
        assert (mol is not None) or (
            smiles is not None
        ), "mol or smiles must be provided"

        self._mol = mol
        self._smiles = smiles
        self._molH = None
        self._is_canon = False

    @property
    def mol(self) -> Type[rdkit.Chem.Mol]:
        if self._mol is None:
            self._mol = rdkit.Chem.MolFromSmiles(self._smiles)
        return self._mol

    @property
    def molH(self) -> Type[rdkit.Chem.Mol]:
        if self._molH is None:
            self._molH = rdkit.Chem.AddHs(self.mol)
        return self._molH

    @property
    def smiles(self) -> str:
        if (self._smiles is None) or not self._is_canon:
            self._smiles = rdkit.Chem.MolToSmiles(self.mol)
        return self._smiles

    @property
    def typePONA(self) -> str:
        if