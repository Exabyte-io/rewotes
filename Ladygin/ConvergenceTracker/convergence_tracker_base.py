import pickle
import os
from driver import Driver


class ConvergenceTracker():
    """Base class for convergence tracker
           Contains methods nessasure for any conv tracker"""
    
    def __init__(self, workdir: str = './', target:str = 'total_energy', encut: float = 600., eps = 1e-2) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)
        """

        self.workdir = workdir
        self.target = target
        self.encut = encut
        self.eps = eps
        self.kpoint_opt = None # optimal kpoint value
        
    def find_opt(self, driver: Driver) -> None:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached"""
        pass

    def save(self, name:str) -> None:
        """save class as self.name.dat"""
        self.name = name
        file = open(self.name+'.dat','wb')
        pickle.dump(self.__dict__, file)
        file.close()

    def load(self) -> None:
        """try load self.name.dat"""

        file = open(self.name+'.dat','rb')
        
        self.__dict__ = pickle.load(file)
    
    def __repr__(self) -> str:
        return f"Convergence Tracker of {self.target}) with eps equal to {self.eps} with results saved at {self.workdir}"


# Usage examples
if __name__ == "__main__":
    Tracker = ConvergenceTracker()
    print(Tracker)
    Tracker.save('TestTracker')
    Tracker.load()
    print(Tracker)
    os.remove('TestTracker.dat')