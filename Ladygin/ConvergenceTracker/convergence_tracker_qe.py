from convergence_tracker_base import ConvergenceTracker
from driver import DriverQE
from utils import consts
from kpoint_scheduler import kpoint_scheduler_uniform

class ConvergenceTrackerQE(ConvergenceTracker):
    """Implimentation of convergence tracker for Quantum Espresso Simulation"""

    def __init__(self, workdir: str = './', target:str = 'total_energy', eps = 1e-2) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)
        """
        super().__init__(workdir, target, encut, eps)
        self.driver = DriverQE(workdir = self.workdir,
                               input_file_name = "pw.in", 
                               encut = 40)


    def _step(self, k_curr: int) -> float:
        """Step of the kpoint scheduler

           k_curr - current k point
        """
        self.driver.gen_input(k_curr)
        self.driver.run()
        return self.driver.extract_target(self.target)
    
    def find_opt(self, driver: DriverQE = None, k_sch: kpoint_scheduler_uniform = kpoint_scheduler_uniform(2, 30)) -> tuple:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
           DriverQE - calculation driver for quantum espresso
        
        """

        if type(driver) != type(None):
            self.driver = driver

        # Ensure working in same directory
        try:
            assert self.driver.workdir == self.workdir:
        except AssertionError:
            self.driver.workdir == self.workdir
        
        self.k_sch = k_sch
        
        # Start with 2 kpoints to estimate error
        target_list = [float('inf')]
        errors = [float('inf')]
        k_curr = float('inf')
        
        while errors[-1] > self.eps * 1e-3 and k_curr > 0:
            k_curr = k_sch.get_next()
            target_curr = self._step(k_curr)
            print(k_curr, target_curr)
            target_list.append(target_curr)
            errors.append(abs(target_list[-1] - target_list[-2]))
            
        if errors[-1] < self.eps:
            self.kpoint_opt = k_curr

        self.target_list = target_list[1:]
        return k_curr, errors[2:]

if __name__ == """__main__""":
    
    k_sch = kpoint_scheduler_uniform(2, 30)
    Tracker = ConvergenceTrackerQE('./', 'total_energy', 40*consts['Ry'])
    
    kpoint_opt, errors = Tracker.find_opt(k_sch = k_sch)
    
    
    print(f"Optimal value of kpoint is {kpoint_opt} with error of {errors[-1]*1e3} meV")