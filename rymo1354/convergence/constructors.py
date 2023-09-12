import numpy as np
from copy import deepcopy
from tabulate import tabulate
from ase.optimize import QuasiNewton
from ase.calculators.vasp import Vasp
from ase.io import read
import sys
import os

class Constructor:
    def __init__(self):
        self._builders = {} 

    def create(self, key, **kwargs):
        builder = self._builders.get(key)
        if not builder:
            raise ValueError(key)
        return builder.build(**kwargs)

class VaspAtoms:
    def __init__(self):
        self._instance = None

    def build(self, directory):
        self._instance = read(os.path.join(directory, 'POSCAR'))
        return self._instance

class EspressoAtoms:
    def __init__(self):
        self._instance = None

    def build(self, directory):
        # Quantum Espresso logic to get atoms
        self._instance = read(os.path.join(directory, '*.in'))
        return self._instance

class AtomsConstructor(Constructor):
    def __init__(self):
        super().__init__()
        self._builders = {'VASP': VaspAtoms(), 
                          'Espresso': EspressoAtoms()} 
    
    def make(self, key, **kwargs):
        return self.create(key, **kwargs)

class CalculatorConstructor(Constructor):
    def __init__(self):
        super().__init__()
        self._builders = {'VASP': VaspCalculator()}

    def make(self, key, **kwargs):
        return self.create(key, **kwargs)

class VaspCalculator:
    def __init__(self):
        self._instance = None

    def build(self, directory, run_command):
        if not self._instance:
            self._instance = Vasp(command=run_command, directory=directory)
            self._instance.read_incar(os.path.join(directory, 'INCAR'))
            self._instance.read_kpoints(os.path.join(directory, 'KPOINTS'))
            self._instance.read_potcar(os.path.join(directory, 'POTCAR'))
        return self._instance

class PropertyConstructor(Constructor):
    def __init__(self):
        super().__init__()
        self._builders = {'Energy': TotalEnergy()}

    def make(self, key, **kwargs):
        return self.create(key, **kwargs)

class TotalEnergy:
    def __init__(self):
        self._instance = None

    def build(self, atoms):
        self._instance = atoms.get_potential_energy()
        return self._instance

class GetFmax(Constructor):
    def __init__(self):
        super().__init__()
        self._builders = {'VASP': VaspFmax()}

    def make(self, key, **kwargs):
        return self.create(key, **kwargs) 

class VaspFmax:
    def __init__(self):
        self._instance = None

    def build(self, calculator):
        ediffg = calculator.asdict().get('inputs', {}).get('ediffg', 1000)
        self._instance = np.abs(ediffg) if np.sign(ediffg) == -1 else ediffg
        return self._instance
    
