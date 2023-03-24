
from abc import ABC

class Calculation(ABC):
    @abstractmethod
    def __init__(self, attribute, **kwargs):
        self.attribute = attribute
        pass 

    @property
    def raw_value(self):
        if self.attribute == 'etotal':
            return self.energy
        elif self.attribute == 'forces':
            return self.forces

class VaspCalculation(Calculation):
    def __init__(self, attribute, **kwargs):
        Calculation.__init__(self, attribute)
        
        from ase.calculators.vasp import Vasp
        if 'path' in kwargs.keys():
            # finished benchmark calculation
            self.atoms = Vasp(restart=True, directory=path).get_atoms()
        elif 'calculator' in kwargs.keys():
            # calculation not yet performed
            self.atoms = calc.get_atoms()

        self.energy = self.atoms.get_potential_energy()
        self.forces = self.atoms.get_forces()
        
    
class QuantumEspressoCalculation(Calculation):
    
