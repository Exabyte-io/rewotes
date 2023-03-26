import numpy as np
import os
from abc import ABC, abstractmethod
from ase.io import read
from error_metric import ErrorMetricScalar, ErrorMetricVector, ErrorMetricMatrix
from utils import PeriodicDftPackages, ConvergenceProperty, ConvergenceParameter
from calculation import VaspCalculation


class ConvergenceTracker(ABC):
    """Holds a bunch of calculations

    Parameters:
    ----------
    parameter: str
    the parameter to converge (e.g. k-points, encut, etc.)

    atoms: ase.atoms.Atoms object
    the material to use

    conv_property: str
    the error metric to use (e.g. total energy, lattice constant, etc.)
    """

    @abstractmethod
    def __init__(
        self,
        path,
        converge_property,
        package=PeriodicDftPackages.VASP,
        cutoff=400,
        eps=1e-2,
        **kwargs,
    ):
        """Constructor; should initialize the series of calculations
        but not start them. Must be overriden
        """

        if converge_property == ConvergenceProperty.etotal:
            self.error_metric = ErrorMetricScalar()
        elif converge_property == ConvergenceProperty.phonon_modes:
            self.error_metric = ErrorMetricVector()

        self.converge_property = converge_property
        self.path = path
        self.package = package
        self.cutoff = cutoff
        self.eps = eps

        # less important parameters
        defined_params = "max_trials"
        for key in defined_params:
            if key in kwargs:
                setattr(self, key, kwargs[key])

        self.read_structure()

    def read_structure(self):
        """Initialize structure."""
        supported_formats = ["POSCAR", "CONTCAR", "poscar", ".vasp", ".cif", ".xyz"]
        files = [
            f
            for f in os.listdir(self.path)
            if os.path.isfile(os.path.join(self.path, f))
        ]
        supported_files = [f for f in files if any([s in f for s in supported_formats])]
        assert (
            len(supported_files) > 0
        ), f"{path} does not contain a supported structure file"
        file_to_read = supported_files[0]
        if len(supported_files):
            print(
                f"Warning: more than one supported files found! Reading the first one: {file_to_read}"
            )
        atoms = read(os.path.join(self.path, file_to_read), index=":")

        if len(atoms) > 1:
            print(
                f"Warning: structure file {file_to_read} has multiple structures, picking the first one."
            )
        atoms = atoms[0]
        self.atoms = atoms

    def read_input(self):
        """initialize other calculation parameters"""
        if self.package == PeriodicDftPackages.vasp:
            self.reference_calculation = VaspCalculation(
                self.converge_property,
                path=self.path,
                atoms=self.atoms,
            )
        elif self.package == PeriodicDftPackages.qe:
            raise NotImplementedError

    @abstractmethod
    def setup_calcs(self):
        pass

    @abstractmethod
    def run_calcs(self, parameter):
        pass


class KpointConvergenceTracker(ConvergenceTracker):
    def __init__(self, path, conv_property, **kwargs):
        ConvergenceTracker.__init__(self, path, conv_property, **kwargs)
        # initialize a series of calculations
        # setting attribute to 'etotal'

    def find_kpoint_series(self):
        cell_lengths = self.atoms.cell.cellpar()[0:3]
        min_ka = 10
        max_ka = 80
        num_ka = 10

        # loop over ka values.
        kpoint_series = []
        for ka in np.linspace(min_ka, max_ka, 10):
            kpoints = []
            for axis_length in cell_lengths:
                kpoints.append(int(ka / axis_length))

            # if already in the set do not add
            if len(kpoint_series) > 0 and kpoints == kpoint_series[-1]:
                continue

            kpoint_series.append(kpoints)
        self.kpoint_series = kpoint_series

    def setup_calcs(self):

        self.find_kpoint_series()

        if self.package == PeriodicDftPackages.vasp:
            calculations = []
            for kpoints in self.kpoint_series:
                calculation = VaspCalculation(
                    self.converge_property,
                    self.atoms,
                    self.path,
                    self.reference_calculation.calculator,
                )
                calculation.kpoints = kpoints
                calculations.append(calculation)

        elif self.package == PeriodicDftPackages.qe:
            raise NotImplementedError

        self.calculations = calculations


class PwCutoffConvergenceTracker(ConvergenceTracker):
    def __init__(self, path, conv_property):
        ConvergenceTracker.__init__(self, path, conv_property)
        raise NotImplementedError
