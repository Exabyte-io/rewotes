"""Base solver class."""
import shutil
from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path


implemented_solvers = ["espresso"]


class BaseSolver(ABC):
    """Base class for solver implementations."""

    def __init__(self, input_dict: dict, input_path: str):
        self.input_dict: dict = deepcopy(input_dict)
        self.solver_path: Path
        self.results_path: Path = Path(input_path)
        self.results: dict = {}
        self.solver_path: Path = self._verify_solver_path()

    @abstractmethod
    def run(self):
        """Run solver-specific executable command."""
        raise NotImplementedError

    @abstractmethod
    def parse_output(self):
        """Implement solver-specific output format parsing."""
        raise NotImplementedError

    def _verify_solver_path(self) -> Path:
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


class SolverNotImplementedError(Exception):
    """Error when the requested solver is not implemented."""


class SolverNotInstalledError(Exception):
    """Error when the requested solver is not installed."""
