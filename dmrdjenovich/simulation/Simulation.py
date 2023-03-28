from computing.executable import Executable
from qe_spec import QESpec
from qe_setup import QESetup
from qe_run import QERun
from qe_process import QEProcess

class Simulation(Executable):
    """
    Class that automatically handles the details of a
    simulation task.
    1) Creates the simulation directory & input files
    2) Runs the simulation executable
    3) Processes the output files
    """
    
    def __init__(self, dir, rsx, setup, run, process):
        self.dir = dir
        self.rsx = rsx
        self.qe_setup = setup
        self.qe_run = run
        self.qe_process = process
        
    @staticmethod
    def construct_default(dir, spec, rsx):
        """
        Creates a default simulation with the given
        specification to run in the given directory
        """
        if isinstance(spec, QESpec):
            return Simulation(dir, rsx, QESetup(dir, spec), QERun(dir, rsx), QEProcess(dir, spec))
        else:
            return None
    
    @staticmethod
    def construct_override(dir, spec, rsx, setup, run, process):
        """
        Creates a simulation with custom setup, run,
        and / or process with the given specification
        to run in the given directory.
        """
        if isinstance(spec, QESpec):
            if setup is None:
                setup = QESetup(dir, spec)
            if run is None:
                run = QERun(dir, rsx)
            if process is None:
                process = QEProcess(dir, spec)
            return Simulation(dir, rsx, setup, run, process)
        else:
            return None
    
    def get_resources(self):
        """
        Returns the minimum computational resources
        required to run this simulation.
        """
        return self.rsx
        
    def run(self, envr):
        """
        Executes this simulation given the computational
        resources passed in.
        """
        if not self.qe_setup.run(envr):
            print("Error detected in setup: " + self.dir)
            return False
        if not self.qe_run.run(envr):
            print("Error detected in run: "  + self.dir)
            return False
        if not self.qe_process.run(envr):
            print("Error detected in process: " + self.dir)
            return False
        return True

