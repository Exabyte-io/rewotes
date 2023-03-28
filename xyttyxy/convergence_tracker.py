import numpy as np
import os
from abc import ABC, abstractmethod
from ase.io import read
from error_metric import ErrorMetricScalar, ErrorMetricVector, ErrorMetricMatrix
from utils import (
    PeriodicDftPackages,
    ConvergenceProperty,
    ConvergenceParameter,
    messages,
)
from calculation import VaspCalculation
from logging import getLogger


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
        workdir,
        path,
        convergence_property,
        package=PeriodicDftPackages.VASP,
        eps=1e-5,
        logger=getLogger("ConvergenceTracker"),
        **kwargs,
    ):
        """Constructor; should initialize the series of calculations
        but not start them. Must be overriden
        """
        if convergence_property == ConvergenceProperty.etotal:
            self.error_metric = ErrorMetricScalar()
        elif convergence_property == ConvergenceProperty.phonon_modes:
            self.error_metric = ErrorMetricVector()

        self.convergence_property = convergence_property
        self.path = path
        self.package = package
        self.eps = eps
        self.logger = logger

        if workdir:
            self.workdir = workdir
        else:
            self.workdir = self.path

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
            self.logger.warning(
                messages("multiple_structure_files").format(file_to_read)
            )

        atoms = read(os.path.join(self.path, file_to_read), index=":")

        if len(atoms) > 1:
            self.logger.warning(
                messages("multiple_structure_images").format(file_to_read)
            )
        atoms = atoms[0]
        self.atoms = atoms.copy()

    def read_input(self):
        """initialize other calculation parameters"""
        if self.package == PeriodicDftPackages.vasp:
            self.reference_calculation = VaspCalculation(
                name="reference",
                conv_property=self.convergence_property,
                path=self.path,
                atoms=self.atoms,
                logger=self.logger,
            )
            # make sure input contains only electronics step related stuff
            # todo: these should be more decoupled from VASP
            self.reference_calculation.no_tetrahedron_smearing()
            self.reference_calculation.set_to_singlepoint()
        elif self.package == PeriodicDftPackages.qe:
            raise NotImplementedError

    @abstractmethod
    def setup_calcs(self):
        pass

    @abstractmethod
    def run_calcs(self):
        error = 1e10
        # first run calculation 0
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        os.chdir(self.workdir)
        calc_low = self.calculations[0]
        calc_low.run()
        for calculation in self.calculations[1:]:
            calc_high = calculation
            calc_high.run()
            error = self.error_metric.error(calc_low, calc_high)
            if error < self.eps:
                # calc_low has the converged mesh
                # since increasing the mesh does not bring significant gain
                self.converged_calc = calc_low
                self.logger.info(messages("converged"))
                os.chdir("..")
                return
            else:
                calc_low = calc_high

        self.logger.warning(messages("unconverged"))
        os.chdir("..")
        # we don't need to record the results separately,
        # they are automatically available in the calculation objects

    @abstractmethod
    def show_results(self):
        pass


class KpointConvergenceTracker(ConvergenceTracker):
    def __init__(self, min_ka, max_ka, num_ka, *args, **kwargs):
        ConvergenceTracker.__init__(self, *args, **kwargs)
        self.min_ka = min_ka
        self.max_ka = max_ka
        self.num_ka = num_ka

    def read_input(self):
        ConvergenceTracker.read_input(self)
        # additionally stop it trying to manipulate kpoints
        self.reference_calculation.no_kspacing()

    def find_kpoint_series(self):
        cell_lengths = self.atoms.cell.cellpar()[0:3]

        # loop over ka values.
        kpoint_series = []
        for ka in np.linspace(self.min_ka, self.max_ka, self.num_ka):
            kpoints = []
            for axis_length in cell_lengths:
                kpoints.append(int(ka / axis_length))

            # if already in the set do not add
            if len(kpoint_series) > 0 and kpoints == kpoint_series[-1]:
                continue

            # if any is 0, do not add
            if any(np.array(kpoints) == 0):
                continue

            kpoint_series.append(kpoints)
        self.kpoint_series = kpoint_series

    def setup_calcs(self):
        self.find_kpoint_series()
        if self.package == PeriodicDftPackages.vasp:
            calculations = []
            for kpoints in self.kpoint_series:
                calculation = VaspCalculation(
                    name="_".join([str(k) for k in kpoints]),
                    conv_property=self.convergence_property,
                    atoms=self.atoms,
                    calculator=self.reference_calculation._calculator,
                    logger=self.logger,
                )
                calculation.kpoints = kpoints
                calculations.append(calculation)

        elif self.package == PeriodicDftPackages.qe:
            raise NotImplementedError

        self.calculations = calculations

    def run_calcs(self):
        ConvergenceTracker.run_calcs(self)
        try:
            self.converged_kpoints = self.converged_calc.kpoints
        except AttributeError:
            self.logger.warning(messages("unconverged_kpts"))
            self.converged_kpoints = [-1, -1, -1]

    def show_results(self):
        xs = []
        ys = []
        labels = []
        for idx, calculation in enumerate(self.calculations):
            if not calculation.calculation_required:
                xs.append(idx)
                ys.append(calculation.etotal)
                label = r"$\times$".join([str(k) for k in calculation.kpoints])
                labels.append(label)
            else:
                break
        import matplotlib.pyplot as plt

        xs = np.array(xs)[0:-1]
        ys = np.abs(np.diff(ys))

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(xs, ys, "bo-")
        ax.set_yscale("symlog", linthresh=self.eps * 1e-2)
        for x, y, label in zip(xs, ys, labels):
            ax.annotate(
                label, (x, y), textcoords="offset points", xytext=(0, 10), ha="center"
            )
        ax.axhline(y=self.eps, color="r", linestyle="dotted")
        ax.annotate(
            "dE", (0, self.eps), textcoords="offset points", xytext=(0, 10), ha="center"
        )

        ax.set_ylabel("Error / eV")
        ax.get_xaxis().set_visible(False)
        plt.savefig(f"{self.workdir}/KpointConvergenceTrack.png")

        self.logger.info(messages("converged_kpts").format(str(self.converged_kpoints)))


class PwCutoffConvergenceTracker(ConvergenceTracker):
    def __init__(self, workdir, path, conv_property, **kwargs):
        ConvergenceTracker.__init__(self, workdir, path, conv_property, **kwargs)
        raise NotImplementedError
