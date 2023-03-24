from abc import ABC

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
    def __init__(self, atoms, conv_property):
        """Constructor; should initialize the series of calculations
        but not start them. Must be overriden
        """
        self.atoms = atoms
        if conv_property == 'etotal':
            self.error_metric = ErrorMetricScalar()
        elif conv_property == 'phonons':
            self.error_metric = ErrorMetricVector()

    def run_calcs(self, parameter):

              
class KpointConvergenceTracker(ConvergenceTracker):
    def __init__(self, atoms, conv_property):
        ConvergenceTracker.__init__(self, atoms, conv_property)

        # initialize a series of calculations,
        # setting attribute to 'etotal'
    
    

class PwCutoffConvergenceTracker(ConvergenceTracker):
    def __init__(self, atoms, conv_property):
        ConvergenceTracker.__init__(self, atoms, conv_property)
        raise NotImplementedError
    
