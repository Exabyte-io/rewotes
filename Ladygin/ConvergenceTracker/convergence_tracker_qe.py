from convergence_tracker_base import ConvergenceTracker
from driver import DriverQE
from utils import consts

class ConvergenceTrackerQE(ConvergenceTracker):
    """Implimentation of convergence tracker for Quantum Espresso Simulation"""

    def __init__(self, workdir: str = './', target:str = 'total_energy', encut: float = 600., eps = 1e-2) -> None:
        """workdir - directory for input files and calculations
           target - target property to optimize
           encut - kinetic energy cutoff (in eV)
           eps - convergence criteria (in meV)
        """
        super().__init__(workdir, target, encut, eps)
    
    def find_opt(self, driver: DriverQE = DriverQE()) -> None:
        """Finds optimal kpoint by running simulations and refining the parameter until convergence reached
           DriverQE - calculation driver for quantum espresso
        
        """
        
        # Start with 2 kpoints to estimate error
        kpoint_list = [2, 4]
        target_list = []
        
        driver.gen_input(self.encut, kpoint_list[0])
        driver.run()
        target_list.append(driver.extract_target(self.target))

        driver.gen_input(self.encut, kpoint_list[1])
        driver.run()
        target_list.append(driver.extract_target(self.target))

        
        error = [abs(target_list[-1] - target_list[-2])]
        
        # while converge go increments of 2
        num_iter = 10
        i = 0
        while error[-1] > self.eps * 1e-3 and i < num_iter:
            i+= 1
            
            kpoint_list.append(kpoint_list[-1] + 2)
            driver.gen_input(self.encut/consts['Ry'], kpoint_list[-1])
            driver.run()
            target_list.append(driver.extract_target(self.target))
            
            error.append(abs(target_list[-1] - target_list[-2]))

        if error[-1] < self.eps:
            self.kpoint_opt = kpoint_list[-1]
        

if __name__ == """__main__""":
    Tracker = ConvergenceTrackerQE('./', 'total_energy', 40*consts['Ry'])

    Tracker.find_opt()

    print(f"Optimal value of kpoint is {Tracker.kpoint_opt}")