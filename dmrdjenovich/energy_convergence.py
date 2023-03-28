from computing.resources import Resources
from convergence import Convergence
import sys

def converge_energy(dir, encut, thresh):
    conv = EnergyConvergence(dir, encut, thresh, Resources(1, -1), homogeneous_k=True)
    if not conv.run(Resources(1, -1)):
        print("Failed")
    return conv.ks[len(conv.ks) - 1]

class EnergyConvergence(Convergence):
    """
    Class responsible for running a series of
    Simulations to do an energy convergence test
    for a given material system.
    """
    
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

def cmd_line():
    print(converge_energy(sys.argv[1], float(sys.argv[2]), float(sys.argv[3])))

cmd_line()
