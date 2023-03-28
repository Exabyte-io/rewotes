from abc import ABC, abstractmethod
from utils import graceful_exit, ConvergenceProperty, messages
from copy import deepcopy as dcopy
import os
import logging


class Calculation(ABC):
    @abstractmethod
    def __init__(self, name, convergence_property, atoms, path, **kwargs):
        self.convergence_property = convergence_property
        self.atoms = atoms.copy()
        self.path = path
        self.name = name

        if "logger" in kwargs.keys():
            self.logger = kwargs["logger"]
        else:
            self.logger = logging.getLogger("Calculation")  # level defaults to warning

    @property
    def raw_value(self):
        if self.convergence_property == ConvergenceProperty.etotal:
            return self.etotal
        elif self.convergence_property == ConvergenceProperty.force:
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
        self, name, convergence_property, atoms, path=None, calculator=None, **kwargs
    ):
        Calculation.__init__(self, name, convergence_property, atoms, path, **kwargs)
        from ase.calculators.vasp import Vasp

        if path:
            _calculator = Vasp()
            try:
                _calculator.read_incar(filename=os.path.join(self.path, "INCAR"))
            except FileNotFoundError:
                self.logger.warning(messages("no_incar").format("INCAR", self.path))
                try:
                    files = [
                        f
                        for f in os.listdir(self.path)
                        if ".incar" in f and os.path.isfile(os.path.join(self.path, f))
                    ]
                    assert isinstance(files, list), "BUG 1"

                    if len(files) == 0:
                        self.logger.warning(
                            messages("no_incar").format(".incar", self.path)
                        )
                        raise FileNotFoundError
                    elif len(files) > 1:
                        self.logger.warning(messages("multiple_incar").format(".incar"))
                        fname = files[0]
                    else:
                        fname = files[0]

                    _calculator.read_incar(filename=os.path.join(self.path, fname))

                except FileNotFoundError:
                    self.logger.error(messages("no_incar_final"))
                    graceful_exit()

        elif calculator:
            _calculator = dcopy(calculator)
        else:
            self.logger.error(messages("no_calculation_input"))
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
            self.logger.warning(messages("set_to_singlepoint"))
        self._calculator.set(ibrion=-1, nsw=0)

    def no_kspacing(self):
        kspacing = self._calculator.float_params["kspacing"]
        if kspacing:
            self.logger.warning(messages("ignore_kspacing"))
        self._calculator.set(kspacing=None)

    @property
    def etotal(self):
        # if self._calculator.calculation_required(self.atoms, ["energy"]):
        if self.calculation_required:
            raise Exception(messages("no_etotal"))
        return self.atoms.get_potential_energy()

    @kpoints.setter
    def kpoints(self, kpts):
        self._calculator.set(kpts=kpts)

    def run(self):
        self.atoms.get_potential_energy()
        self.calculation_required = False
