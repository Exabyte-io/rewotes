

class ConvTracker:
    """
    Args:
        kwargs:
            cutoff (float): Desired energy cutoff in eV.
            energy (list): Total energy values. Can be used as a pseudo-restart to convergence.

    Attributes:
        cutoff (float): Desired energy cutoff in eV.
        energy (list): List of energy values used to check for convergence.
    """

    def __init__(self, config, job_endpoints, cutoff=1e-5, energy=[]):
        self.cutoff = cutoff            # Units = eV
        self.energy = energy            # Array of energies can be passed in to continue a job set.

    def check_convergence(self):
        """
        Check if energy convergence reached.
        """
        if len(self.energy) < 2:
            return False
        else:
            return abs(self.energy[-1] - self.energy[-2]) <= self.cutoff

