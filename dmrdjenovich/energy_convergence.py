from convergence import Convergence

class EnergyConvergence(Convergence):
    """
    Class responsible for running a series of
    Simulations to do a energy convergence test
    for a given material system.
    """
    
    INPUT_NAME = "input.txt"
    
    def __init__(self, dir, encut, thresh, rsx, input_name="input.txt", meta=None, obtusify=True, homogeneous_k=False):
        """
        Sets up a convergence test in the specified
        directory, using the provided QESpec, and
        a metadata string for logging purposes.
        """
        super(EnergyConvergence, self).__init__(dir, encut, thresh, rsx, input_name=input_name, meta=meta, obtusify=obtusify, homogeneous_k=homogeneous_k)
      
    def get_results(self, analysis):
        energies = analysis.get_energy()
        return energies[len(energies) - 1]
