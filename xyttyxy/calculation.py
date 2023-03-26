from abc import ABC, abstractmethod
from utils import graceful_exit


class Calculation(ABC):
    @abstractmethod
    def __init__(self, conv_property, atoms, path, **kwargs):
        self.conv_property = conv_property
        self.atoms = atoms
        self.path = path

    @property
    def raw_value(self):
        if self.conv_property == ConvergenceProperty.etotal:
            return self.energy
        elif self.conv_property == ConvergenceProperty.force:
            return self.forces

    @property
    @abstractmethod
    def kpoints(self):
        pass

    @kpoints.setter
    @abstractmethod
    def kpoints(self, kpts):
        pass


class VaspCalculation(Calculation):
    def __init__(self, conv_property, atoms, path=None, calculator=None, **kwargs):
        Calculation.__init__(self, conv_property, atoms, path, **kwargs)

        from ase.calculators.vasp import Vasp

        if path:
            calculator = Vasp()
            try:
                calculator.read_incar(filename=f"{self.path}/INCAR")
                self.calculator = calculator
            except FileNotFoundError:
                print(f"No INCAR exists at path {self.path} specified")
                graceful_exit()

        elif calculator:
            self.calculator = calculator

    @property
    def kpoints(self):
        return self.calculator.kpts

    @kpoints.setter
    def kpoints(self, kpts):
        self.calculator.set(kpts=kpts)
