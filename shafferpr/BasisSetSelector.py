from ase import Atoms
from ase.calculators import nwchem
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from rdkit import Chem
from rdkit.Chem import AllChem


class Molecule(object):
    def __init__(self,input_file_path='',smiles_string=''):
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

    def calculate_energy(self,basis_set):
        calc=nwchem.NWChem(
              dft=dict(maxiter=2000,
                       xc='B3LYP'),
              basis=basis_set)
        self.atoms.calc=calc
        #BFGS(self.atoms).run(fmax=0.01)
        energy=self.atoms.get_potential_energy()
        return energy
    
    def calculate_vibrational_frequency(self,basis_set):
        calc=nwchem.NWChem(
              dft=dict(maxiter=2000,
                       xc='B3LYP'),
              basis=basis_set)
        self.atoms.calc=calc
        BFGS(self.atoms).run(fmax=0.01)
        vib = Vibrations(self.atoms)
        vib.run()
        vibrational_energy=vib.get_energies()[0]
        return vibrational_energy

    






calc = nwchem.NWChem(label='calc/nwchem',
              dft=dict(maxiter=2000,
                       xc='B3LYP'),
              basis='6-31+G*')