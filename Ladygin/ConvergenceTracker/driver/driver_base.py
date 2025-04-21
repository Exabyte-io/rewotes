from .calculator import calculator


class Driver():
    """Base class for running simulation"""

    def __init__(self, workdir: str, calculator:calculator, input_file_name: str, encut: float) -> None:
        """input_file - simulation settings input file name
           workdir - working directory
           calculator - job runner for point calculation
           encut - energy cutoff for point calc
        """
        pass


    def gen_input(self, kpoint: float):
        """Generates the input file for the driver based on input parameters
        
            kpoint - kpoint dimension    
        """
        pass

    def calc(self) -> None:
        """Runs point simulation in the current folder"""
        pass


    def extract_target(self, target:str) -> float:
        """Extracts the target from the simulation output"""
        pass
