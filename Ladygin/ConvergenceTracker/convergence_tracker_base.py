import pickle
import os
from .driver import Driver
from .search import kpoint_scheduler

class ConvergenceTracker():
    """Base class for convergence tracker
           Contains methods nessasure for any conv tracker"""
    
    def __init__(self, workdir: str, target:str, eps: float, driver: Driver) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)
        """

        self.workdir = workdir
        self.target = target
        self.eps = eps
        self.kpoint_opt = None # optimal kpoint value

    def save_stat(self, k_list: list, target_list: list, errors: list) -> None:
        """Save stats of kpoint convergence iterations"""
        pass
        
    def _step(self, k_curr: int) -> float:
        """Step of the kpoint scheduler

           k_curr - current k point
        """
        pass
        
    def find_opt(self, driver: Driver, k_sch: kpoint_scheduler) -> None:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
        
           driver - drives for simulations
           k_sch - kpoint generator
        """
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