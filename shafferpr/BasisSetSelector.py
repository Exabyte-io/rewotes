from ase import Atoms
from ase.calculators import nwchem
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os

class Molecule(object):
    def __init__(self,input_file_path='',smiles_string='', label='', reference_value=None):
        self.label=label
        if input_file_path is not '':
            self.atoms=self.atoms_from_file(input_file_path)
        elif smiles_string is not '':
            self.smiles_string=smiles_string
            self.atoms=self.atoms_from_smiles_string(smiles_string)
        else:
            raise RuntimeError("You must specify either an input file or a smiles string to creat a Molecule")

    def atoms_from_smiles_string(self,smiles_string):
        mol=Chem.MolFromSmiles(smiles_string)
        AllChem.Compute2DCoords(mol)
        conf=mol.GetConformer()
        positions=conf.GetPositions()
        atoms=mol.GetAtoms()
        atomic_symbols=[x.GetSymbol() for x in atoms]
        atoms=Atoms(symbols=atomic_symbols,positions=positions)
        return atoms

    def atoms_from_file(self,input_file_path):
        with open(input_file_path,"r") as f:
            lines=f.readlines()[2:]
        atomic_symbols=[x.split()[0] for x in lines]
        positions=[[float(y) for y in x.split()[1:]] for x in lines]
        atoms=Atoms(symbols=atomic_symbols,positions=positions)
        return atoms


    

class Calculator(object):
    def __init__(self,basis_set='STO-2G',property='energy',use_reference_calculator=False):
        if use_reference_calculator:
             self.calc=nwchem.NWChem(
                theory='ccsd'
            )
        else:
            self.calc=nwchem.NWChem(
                dft=dict(maxiter=2000,
                        xc='B3LYP'
                ),
                basis=basis_set
            )

        self.property=property
        if self.property == 'energy':
            self.calculate_property=self.calculate_energy

        elif self.property == 'vibrational_energy':
            self.calculate_property=self.calculate_vibrational_energy
    
    def optimize(self,molecule):
        molecule.atoms.calc=self.calc
        BFGS(molecule.atoms).run(fmax=0.01)

    def calculate_energy(self,molecule):
        molecule.atoms.calc=self.calc
        energy=molecule.atoms.get_potential_energy()
        return energy

    def calculate_vibrational_energy(self,molecule):
        molecule.atoms.calc=self.calc
        BFGS(molecule.atoms).run(fmax=0.01)
        vib = Vibrations(molecule.atoms)
        vib.run()
        vibrational_energy=vib.get_energies()

        os.system("rm vib.*")
        return np.abs(vibrational_energy[0])


