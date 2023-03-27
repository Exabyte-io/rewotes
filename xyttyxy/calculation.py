from abc import ABC, abstractmethod
from utils import graceful_exit, ConvergenceProperty
from copy import deepcopy as dcopy


class Calculation(ABC):
    @abstractmethod
    def __init__(self, name, conv_property, atoms, path, **kwargs):
        self.conv_property = conv_property
        self.atoms = atoms.copy()
        self.path = path
        self.name = name

    @property
    def raw_value(self):
        if self.conv_property == ConvergenceProperty.etotal:
            return self.etotal
        elif self.conv_property == ConvergenceProperty.force:
            raise NotImplementedError

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
    def __init__(
        self, name, conv_property, atoms, path=None, calculator=None, **kwargs
    ):
        Calculation.__init__(self, name, conv_property, atoms, path, **kwargs)

        from ase.calculators.vasp import Vasp

        if path:
            _calculator = Vasp()
            try:
                _calculator.read_incar(filename=f"{self.path}/INCAR")

            except FileNotFoundError:
                print(f"No INCAR exists at path {self.path} specified")
                graceful_exit()
        elif calculator:
            _calculator = dcopy(calculator)
        else:
            print("Either path to an INCAR or a _calculator object must be supplied")
            graceful_exit()

        _calculator.set(directory=name)
        self._calculator = _calculator
        self.atoms.calc = _calculator
        self.calculation_required = True

    @property
    def kpoints(self):
        return self._calculator.kpts

    def no_tetrahedron_smearing(self):
        ismear = self._calculator.int_params["ismear"]
        if ismear == -5:
            self._calculator.set(ismear=0, sigma=0.02)

    def set_to_singlepoint(self):
        ibrion = self._calculator.int_params["ibrion"]
        nsw = self._calculator.int_params["nsw"]
        if ibrion != -1 and nsw != 0:
            print("setting to a single point calculation for convergence tracking")
        self._calculator.set(ibrion=-1, nsw=0)

    @property
    def etotal(self):
        # this should be tested as a read-only attribute
        if self._calculator.calculation_required(self.atoms, ["energy"]):
            print("Total energy not available yet, call VaspCalculation.run first")
        return self.atoms.get_potential_energy()

    @kpoints.setter
    def kpoints(self, kpts):
        self._calculator.set(kpts=kpts)

    def run(self):
        self.atoms.get_potential_energy()
        self.calculation_required = False
