from abc import ABC, abstractmethod


from .calculator import calculator


class Driver(ABC):
    """Base class for running simulation"""


    @abstractmethod
    def __init__(self, workdir: str, calculator:calculator, input_file_name: str, encut: float, **kwargs) -> None:
        """input_file_name - simulation settings input file name
           workdir - working directory
           calculator - job runner for point calculation
           encut - energy cutoff for point calc
        """
        self.workdir = workdir
        self.calculator = calculator
        self.input_file_name = input_file_name
        self.encut =  encut

    @abstractmethod
    def gen_input(self, kpoint: float):
        """Generates the input file for the driver based on input parameters
        
            kpoint - kpoint dimension    
        """
        raise NotImplementedError

    @abstractmethod
    def calc(self) -> None:
        """Runs point simulation in the current folder"""
        raise NotImplementedError

    @abstractmethod
    def extract_target(self, target:str) -> float:
        """Extracts the target from the simulation output"""
        raise NotImplementedError
