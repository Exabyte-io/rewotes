from abc import ABC, abstractmethod
from error_metric import ErrorMetricScalar, ErrorMetricVector, ErrorMetricMatrix
from utils import PeriodicDftPackages, ConvergenceProperty, ConvergenceParameter

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
    def __init__(self, atoms, conv_property, package = PeriodicDftPackages.VASP, cutoff = 400, eps = 1e-2, **kwargs):
        """Constructor; should initialize the series of calculations
        but not start them. Must be overriden
        """
        self.atoms = atoms
        if conv_property == ConvergenceProperty.total_energy:
            self.error_metric = ErrorMetricScalar()
        elif conv_property == ConvergenceProperty.phonon_modes:
            self.error_metric = ErrorMetricVector()
            
        self.package = package
        self.cutoff = cutoff
        self.eps = eps

        # less important parameters
        defined_params = ('max_trials')
        for key in defined_params:
            if key in kwargs:
                setattr(self, key, kwargs[key])

    def setup_calcs(self):
        if self.package == PeriodicDftPackages.VASP:
            pass

    def run_calcs(self, parameter):
        pass


class KpointConvergenceTracker(ConvergenceTracker):
    def __init__(self, atoms, conv_property, **kwargs):
        ConvergenceTracker.__init__(self, atoms, conv_property, **kwargs)
        # initialize a series of calculations
        # setting attribute to 'etotal'


class PwCutoffConvergenceTracker(ConvergenceTracker):
    def __init__(self, atoms, conv_property):
        ConvergenceTracker.__init__(self, atoms, conv_property)
        raise NotImplementedError
