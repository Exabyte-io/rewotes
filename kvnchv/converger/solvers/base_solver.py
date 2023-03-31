"""Base solver class."""
import shutil
from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path

from ..exceptions import (SolverNotInstalledError,
                          UnsupportedParameterError)

# solver name and executable name on path
implemented_solvers = {"espresso": "pw.x"}


class BaseSolver(ABC):
    """Base class for solver implementations.

    Attributes
    ----------
    input_dict : dict
        Contains name and solver_path
    input_path : str
        Directory containing expected input structure
    parameter_set : dict
        Parameter name(s) and value(s) to replace in the provided input file
    supported_parameters : list
        Valid parameterizable inputs, correspond to keys in parameter_set
    results : dict
        All implemented return values from simulation. Keys match 'target'
    solver_path : Path
        Resolved solver path, provided or default (on PATH)
    """

    def __init__(self, input_dict: dict, input_path: str, parameter_set: dict):
        """Construct base Solver object.

        Parameters
        ----------
        input_dict : dict
            Contains name and solver_path
        input_path : str
            Directory containing expected input structure
        parameter_set : dict
            Parameter name(s) and value(s) to replace in the provided input file
        """
        self.input_dict: dict = deepcopy(input_dict)
        self.input_path: Path = input_path
        self.parameter_set: dict = parameter_set
        self.supported_parameters: list = []
        self.results: dict = {}
        self.solver_path: Path = self._validate_solver_path()

    @abstractmethod
    def run(self):
        """Run solver-specific executable command."""
        raise NotImplementedError

    @abstractmethod
    def parse_output(self, outdir: Path):
        """Implement solver-specific output format parsing."""
        raise NotImplementedError

    def _validate_solver_path(self) -> Path:
        """Check solver is installed locally."""
        if "solver_path" in self.input_dict:
            try_path = self.input_dict["solver_path"]
            # fallback to check for solver on PATH
            if not Path(try_path).exists():
                try_path = shutil.which(implemented_solvers[self.input_dict["name"]])
        else:
            try_path = shutil.which(implemented_solvers[self.input_dict["name"]])

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
