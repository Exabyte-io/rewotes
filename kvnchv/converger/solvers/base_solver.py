"""Base solver class."""
import shutil
from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path


implemented_solvers = ["espresso"]


class BaseSolver(ABC):
    """Base class for solver implementations."""

    def __init__(self, input_dict: dict, input_path: str, parameter_set: dict):
        self.input_dict: dict = deepcopy(input_dict)
        self.input_path: Path = input_path
        self.parameter_set: dict = parameter_set
        self.solver_path: Path
        self.supported_parameters: list
        self.results_path: Path = Path(input_path)
        self.results: dict = {}
        self.solver_path: Path = self._validate_solver_path()

    @abstractmethod
    def run(self):
        """Run solver-specific executable command."""
        raise NotImplementedError

    @abstractmethod
    def parse_output(self):
        """Implement solver-specific output format parsing."""
        raise NotImplementedError

    def _validate_solver_path(self) -> Path:
        """Check solver is implemented in package and installed locally."""
        if self.input_dict["name"] not in implemented_solvers:
            raise SolverNotImplementedError

        if "solver_path" in self.input_dict:
            try_path = self.input_dict["solver_path"]
        else:
            try_path = shutil.which("pw.x")
            if not try_path:
                raise SolverNotInstalledError

        if not Path(try_path).exists():
            raise SolverNotInstalledError
        return Path(try_path)

    def _validate_parameters(self):
        """Check editable parameters are supported."""
        if not all(param in self.supported_parameters for param in self.parameter_set):
            raise UnsupportedParameterError

    @abstractmethod
    def _generate_new_file(self, new_input_fname: str):
        """Implement solver-specific input file parameter replacement."""
        raise NotImplementedError


class SolverNotImplementedError(Exception):
    """Error when the requested solver is not implemented."""


class SolverNotInstalledError(Exception):
    """Error when the requested solver is not installed."""


class UnsupportedParameterError(Exception):
    """Error when provided parameters are unsupported."""
