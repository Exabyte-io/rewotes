from .convergence_tracker_base import ConvergenceTracker
from .driver import DriverQE
from .utils import consts
from .kpoint_scheduler import kpoint_scheduler_uniform

import os
import warnings
from .exceptions import NotConvergedWarning


class ConvergenceTrackerQE(ConvergenceTracker):
    """Implimentation of convergence tracker for Quantum Espresso Simulation"""

    def __init__(self, workdir: str = './', target:str = 'total_energy', eps:float = 1e-2) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)
        """
        self.workdir = workdir
        self.target = target
        self.eps = eps
        self.driver = DriverQE(workdir = self.workdir,
                               input_file_name = "pw.in", 
                               encut = 40)



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
    
    def find_opt(self, driver: DriverQE = None, k_sch: kpoint_scheduler_uniform = kpoint_scheduler_uniform(2, 30)) -> tuple:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
           DriverQE - calculation driver for quantum espresso

           Returns:
           k_list - list of kpoints
           target_list - list of target properties
           errors - list of error in target properties
        """

        if type(driver) != type(None):
            self.driver = driver

        # Ensure working in same directory
        try:
            assert self.driver.workdir == self.workdir
        except AssertionError:
            self.driver.workdir = self.workdir
        
        self.k_sch = k_sch
        
        # Start with 2 kpoints to estimate error
        target_list = [float('inf')]
        errors = [float('inf')]
        k_curr = float('inf')
        k_list = []

        it = 1
        print(f"iter, k_curr, {self.target}, error")
        while errors[-1] > self.eps * 1e-3 and k_curr > 0:
            k_curr = k_sch.get_next()
            k_list.append(k_curr)
            
            target_curr = self._step(k_curr)
            
            target_list.append(target_curr)
            errors.append(abs(target_list[-1] - target_list[-2]))
            
            print("%d %d %.6f %.6f" % (it, k_curr, target_curr, errors[-1]))
            it += 1
            
        if errors[-1] > self.eps:
            warnings.warn("The search procedure not converged. Result will return last kpoint and error", NotConvergedWarning)

        # Saving statistics
        self.save_stat(k_list[1:], target_list[2:], errors[2:])
        
        return k_list[1:], target_list[2:], errors[2:]

if __name__ == """__main__""":
    
    k_sch = kpoint_scheduler_uniform(2, 30)
    Tracker = ConvergenceTrackerQE('./', 'total_energy', 40*consts['Ry'])
    
    k_list, target_list, errors = Tracker.find_opt(k_sch = k_sch)
    
    
    print(f"Optimal/Final value of kpoint is {k_list[-1]} with error of {errors[-1]*1e3} meV")