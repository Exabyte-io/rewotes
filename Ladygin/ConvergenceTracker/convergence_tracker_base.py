from abc import ABC, abstractmethod
import pickle
import os


from .driver import Driver
from .search import kpoint_scheduler, metrix



class ConvergenceTracker(ABC):
    """Base class for convergence tracker
           Contains methods nessasure for any conv tracker"""

    @abstractmethod
    def __init__(self, workdir: str,
                 target:str,
                 eps: float,
                 driver: Driver,
                 k_sch: kpoint_scheduler,
                 metrix: metrix,
                 **kwargs) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)

           driver - driver for calculations
           k_sch - kpoint scheduler for seach kpoint gen
           metrix - metrix to compute an error
        """

        self.workdir = workdir
        self.target = target
        self.eps = eps
        self.kpoint_opt = None # optimal kpoint value

        self.driver = driver
        self.k_sch = k_sch
        self.metrix = metrix

    @abstractmethod
    def save_stat(self, k_list: list, target_list: list, errors: list) -> None:
        """Save stats of kpoint convergence iterations"""
        raise NotImplementedError

    @abstractmethod
    def _step(self, k_curr: int) -> float:
        """Step of the kpoint scheduler

           k_curr - current k point
        """
        raise NotImplementedError

    @abstractmethod
    def find_opt(self) -> None:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
        
           driver - drives for simulations
           k_sch - kpoint generator
        """
        raise NotImplementedError

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
