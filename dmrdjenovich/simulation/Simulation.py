from computing.executable import Executable
from simulation.qe_spec import QESpec

class Simulation(Executable):
    """
    Class that automatically handles the details of a
    simulation task.
    1) Creates the simulation directory & input files
    2) Runs the simulation executable
    3) Processes the output files
    """
    
    def __init__(self, rsx, dir, setup, run, process):
        self.dir = dir
        self.rsx = rsx
        self.setup = setup
        self.run = run
        self.process = process
        
    @staticmethod
    def construct_default(dir, spec):
        """
        Creates a default simulation with the given
        specification to run in the given directory
        """
        if isinstance(spec, QESpec):
            return Simulation(dir, QESetup(dir, spec), QERun(dir, spec), QEProcess(dir))
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
            return Simulation(setup, run, process)
        else:
            return None
    
    def get_resources(self):
        """
        Returns the minimum computational resources
        required to run this simulation.
        """
        return self.rsx
        
    def exec(self, envr):
        """
        Executes this simulation given the computational
        resources passed in.
        """
        if not self.setup.exec(envr):
            print("Error detected in setup: " + self.dir)
            return False
        if not self.run.exec(envr):
            print("Error detected in run: "  + self.dir)
            return False
        if not self.process.exec(envr):
            print("Error detected in process: " + self.dir)
            return False
        return True
