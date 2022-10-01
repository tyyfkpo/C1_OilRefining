import pandas as pd

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors,Descriptors,Draw,Fragments
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import MACCSkeys

from itertools import groupby

import sys


import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 300,300

from typing import Dict, Iterator, Type
try:
    import flask
    from flask import Markup
except ImportError:
    flask = None
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from collections import defaultdict




class Molecule:
    all_molecules = []

    
    def __init__(self, mol: Type[rdkit.Chem.Mol] = None, smiles: str = None) -> None:
        assert (mol is not None) or (
            smiles is not None
        ), "mol or smiles must be provided"

        self.__mol = mol
        self.__smiles = smiles
        self.__molH = None
        self.__is_canon = False
        Molecule.all_molecules.append(self.smiles)

    @property
    def mol(self) -> Type[rdkit.Chem.Mol]:
        if self.__mol is None:
            self.__mol = rdkit.Chem.MolFromSmiles(self.__smiles)
        return self.__mol

    @property
    def molH(self) -> Type[rdkit.Chem.Mol]:
        if self.__molH is None:
            self.__molH = rdkit.Chem.AddHs(self.mol)
        return self.__molH

    @property
    def smiles(self) -> str:
        if (self.__smiles is None) or not self.__is_canon:
            self.__smiles = rdkit.Chem.MolToSmiles(self.mol)
        return self.__smiles

    def generate_atom_bond(self) -> list:
        __mol_bond_info = []
        for i, bond in enumerate(self.mol.GetBonds()):
            # print(i,bond)
            __mol_bond_info.append(bond.GetBondType().name)
        return __mol_bond_info

    @property
    def typePONA(self) -> str:
        if Chem.rdMolDescriptors.CalcNumAromaticRings(self.mol):
            self.__typePONA = 'A'
            
        elif Chem.rdMolDescriptors.CalcNumAliphaticRings(self.mol):
            self.__typePONA = 'N'
            
        elif 'DOUBLE' in self.generate_atom_bond():
            self.__typePONA = 'O'
            
        else:
            self.__typePONA = 'P'
            
        return self.__typePONA
    
    @property
    def numANrings(self) -> dict:
        __numANrings = {'Arings':0, 'Nrings':0}
        if (self.typePONA == 'A') or  (self.typePONA == 'N'):
            numArings = Chem.rdMolDescriptors.CalcNumAromaticRings(self.mol)
            numNrings = Chem.rdMolDescriptors.CalcNumAliphaticRings(self.mol)
            __numANrings = {'Arings':numArings, 'Nrings':numNrings}       
        return __numANrings
    
    @property
    def mw(self) -> float:
        self._mw = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
        return self._mw
    
    @property
    def formula(self) -> str:
        self._formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
        return self._formula
    
    def ele_comp(self) -> dict:
        ele_comp_dict =  defaultdict(lambda : 0)
        for atom in self.molH.GetAtoms():
            ele_comp_dict[atom.GetSymbol()] += 1  # GetSymbol()
        return dict(ele_comp_dict)
    
    def save_image(self, label = 'atomNote', showindices = True, figsize=(300, 300), filename = 'temp.svg'):
        for atom in self.mol.GetAtoms():
            atom.SetProp(label, str(atom.GetIdx()+1))
        #Draw.MolToFile(self.mol, filename = 'temptest.png')            
            
        rdDepictor.Compute2DCoords(self.mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
        if not showindices:
            drawer.drawOptions().addAtomIndices = False
        drawer.SetFontSize(0.8)
        drawer.DrawMolecule(self.mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
    
        if flask:
            return Markup(svg)
        else:
            with open(filename, "w") as file:
                file.write(svg)    
            #return svg
    
    def smf_A(self) -> dict:
        # eignvalues need to be retained, in case only det can't differentiate some fetaures.
        pass
    
    
        
'''  
smiles_OX = 'Cc1ccccc1C'
mol_OX = Molecule(smiles = smiles_OX)
mol_OX.mol
mol_OX.smiles
mol_OX.molH
mol_OX.typePONA
mol_OX.ele_comp()

smiles_O = 'C=CCC(C)CC'
mol_O = Molecule(smiles=smiles_O)
mol_O.typePONA
mol_O.mol
mol_O.smiles
mol_O.molH
'''

# =============================================================================
# class AromaticMol(Molecule):
#     def __init__(self, mol = None, smiles = None):
#         super().__init__(mol, smiles)
#         # self.__numArings = None
#     @property
#     def numArings(self):
#         self.__numArings = Chem.rdMolDescriptors.CalcNumAromaticRings(self.mol)
#         return self.__numArings
# 
# 
# smiles_A = 'C1=CC=C2C=C(C=CC2=C1)C(=S)N'
# mol_A = AromaticMol(smiles=smiles_A)
# mol_A.smiles
# mol_A.typePONA
# mol_A.numANrings
# 
# smiles_A2 = 'CCC(C)CCCC1=CC=C2C=C(C=CC2=C1)C(=S)N'
# mol_A2 = AromaticMol(smiles=smiles_A2)
# mol_A2.smiles
# mol_A2.typePONA
# mol_A2.numANrings
# 
# 
# 
# def smiles_draw(smiles_list):
#     mols_list = []
# 
#     for smi in smiles_list:
#         mol = Chem.MolFromSmiles(smi)
#         mols_list.append(mol)
# 
#     img = Draw.MolsToGridImage(
#         mols_list,
#         molsPerRow=4,
#         subImgSize=(200, 200),
#         legends=['' for x in mols_list])
#     return img
# =============================================================================




