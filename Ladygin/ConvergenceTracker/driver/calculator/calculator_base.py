from abc import ABC, abstractmethod

class calculator(ABC):
    """Class for point calculations"""

    @abstractmethod
    def __init__(self, path_to_exec: str, **kwargs) -> None:
        """path_to_exec - path to executable"""
        self.path_to_exec = path_to_exec

    @abstractmethod
    def run_scf(self, input_file:str, output_file:str) -> None:
        """Runs scf calculation"""
        raise NotImplementedError
