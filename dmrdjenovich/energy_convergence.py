import convergence

class EnergyConvergence(Convergence):
    """
    Class responsible for running a series of
    Simulations to do a energy convergence test
    for a given material system.
    """
    
    INPUT_NAME = "input.txt"
    
    def __init__(self, dir, encut, rsx, meta=None, obtusify=True, homogeneous_k=False):
        """
        Sets up a convergence test in the specified
        directory, using the provided QESpec, and
        a metadata string for logging purposes.
        """
        super(Convergence, self).__init__(dir, encut, rsx, meta=meta, obtusify=obtusify, homogeneous_k=homogeneous_k)
      
    @abstractmethod
    def get_results(self, analysis):
        energies = analysis.get_energy()
        return energies[len(energies) - 1]
