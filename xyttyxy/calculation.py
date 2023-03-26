from abc import ABC, abstractmethod
from utils import graceful_exit


class Calculation(ABC):
    @abstractmethod
    def __init__(self, name, conv_property, atoms, path, **kwargs):
        self.conv_property = conv_property
        self.atoms = atoms
        self.path = path
        self.name = name

    @property
    def raw_value(self):
        if self.conv_property == ConvergenceProperty.etotal:
            return self.energy
        elif self.conv_property == ConvergenceProperty.force:
            return self.forces

    @property
    @abstractmethod
    def etotal(self):
        pass

    @property
    @abstractmethod
    def kpoints(self):
        pass

    @kpoints.setter
    @abstractmethod
    def kpoints(self, kpts):
        pass

    @abstractmethod
    def run(self):
        pass

class VaspCalculation(Calculation):
    def __init__(self, name, conv_property, atoms, path=None, calculator=None, **kwargs):
        Calculation.__init__(self, name, conv_property, atoms, path, **kwargs)

        from ase.calculators.vasp import Vasp

        if path:
            calculator = Vasp()
            try:
                calculator.read_incar(filename=f"{self.path}/INCAR")
            except FileNotFoundError:
                print(f"No INCAR exists at path {self.path} specified")
                graceful_exit()

        calculator.set(directory = name)
        self.calculator = calculator
        self.atoms.calc = calculator
        self.calculation_required = True

    @property
    def kpoints(self):
        return self.calculator.kpts

    @property
    def etotal(self):
        # this should be tested as a read-only attribute
        if self.calculator.calculation_required(self.atoms, ['energy']):
            print('Total energy not available yet, call VaspCalculation.run first')
        return self.atoms.get_potential_energy()

    @kpoints.setter
    def kpoints(self, kpts):
        self.calculator.set(kpts=kpts)

    def run(self):
        self.atoms.get_potential_energy()
        self.calculation_required = False

    
