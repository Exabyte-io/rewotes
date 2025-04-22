from .convergence_tracker_base import ConvergenceTracker

import os
import warnings

from .driver import DriverQE
from .search import kpoint_scheduler, kpoint_scheduler_uniform
from .search import metrix, mae

from .utils.exceptions import NotConvergedWarning
from .utils.consts import consts

class ConvergenceTrackerQE(ConvergenceTracker):
    """Implimentation of convergence tracker for Quantum Espresso Simulation"""

    def __init__(self, workdir: str = './', target:str = 'total_energy', eps:float = 1e-2, 
                 driver: DriverQE = DriverQE(),
                 k_sch: kpoint_scheduler = kpoint_scheduler_uniform(2, 30),
                 metrix: metrix = mae(),
                 **kwargs) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)

           driver - driver for calculations
           k_sch - kpoint scheduler for seach kpoint gen
           metrix - metrix to compute an error
        """
        super().__init__(workdir, target, eps, driver, k_sch, metrix, **kwargs)


    def save_stat(self, k_list: list, target_list: list, errors: list) -> None:
        """Save stats of kpoint convergence iterations"""

        out_path = os.path.join(self.workdir, 'stats.dat')
        with open(out_path, 'w') as f:
            f.write(f"kpoint {self.target} error\n")
            for i in range(len(k_list)):
                f.write('%d %.8f %.8f\n' % (k_list[i], target_list[i], errors[i]))
        
    def _step(self, k_curr: int) -> float:
        """Step of the kpoint scheduler

           k_curr - current k point
        """
        self.driver.gen_input(k_curr)
        self.driver.calc()
        return self.driver.extract_target(self.target)
    
    def find_opt(self) -> tuple:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
            
           Returns:
           k_list - list of kpoints
           target_list - list of target properties
           errors - list of error in target properties
        """

        # Ensure working in same directory
        try:
            assert self.driver.workdir == self.workdir
        except AssertionError:
            self.driver.workdir = self.workdir
        
        # Start with 2 kpoints to estimate error
        target_list = [float('inf')]
        errors = [float('inf')]
        k_curr = float('inf')
        k_list = []

        it = 1
        print(f"iter, k_curr, {self.target}, error")
        while errors[-1] > self.eps * 1e-3:
            k_curr = self.k_sch.get_next()
            if k_curr < 0:
                break
            k_list.append(k_curr)
            
            target_curr = self._step(k_curr)
            
            target_list.append(target_curr)
            errors.append(self.metrix(target_list[-1], target_list[-2]))
            
            print("%d %d %.6f %.6f" % (it, k_curr, target_curr, errors[-1]))
            it += 1
            
        if errors[-1] > self.eps * 1e-3:
            warnings.warn("The search procedure not converged. Result will return last kpoint and error", NotConvergedWarning)

        # Saving statistics
        self.save_stat(k_list[1:], target_list[2:], errors[2:])
        
        return k_list[1:], target_list[2:], errors[2:]