from abc import ABC
from utils import graceful_exit


class Calculation(ABC):
    @abstractmethod
    def __init__(self, conv_property, path, **kwargs):
        self.conv_property = conv_property
        self.path = path

    @property
    def raw_value(self):
        if self.conv_property == ConvergenceProperty.etotal:
            return self.energy
        elif self.conv_property == ConvergenceProperty.force:
            return self.forces


class VaspCalculation(Calculation):
    def __init__(self, conv_property, path=None, calculator=None, **kwargs):
        Calculation.__init__(self, conv_property, path, **kwargs)

        from ase.calculators.vasp import Vasp

        if path:
            calc = Vasp()
            try:
                calc.read_incar(filename=f"{self.path}/INCAR")
            except FileNotFoundError:
                print(f"No INCAR exists at path {self.path} specified")
                graceful_exit()

        elif calculator:
            calc = calculator

        self.energy = self.atoms.get_potential_energy()
        self.forces = self.atoms.get_forces()
