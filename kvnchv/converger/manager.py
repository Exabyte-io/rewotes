"""Convergence workflow manager."""

from copy import deepcopy
from pathlib import Path

from jsonschema import ValidationError, validate

from .schema import input_schema
from .solvers import EspressoSolver

solvers = {"espresso": EspressoSolver}


class Manager(object):
    """Main class for convergence workflow run."""

    def __init__(self, input_dict: dict):
        self.input_dict: dict = deepcopy(input_dict)
        self.schema: dict = input_schema

        self.solver_input: dict
        self.input_path: Path
        self.tol: float
        self.target: str

        self.results: list = []

        self.validate_input()
        self._parse_input()
        self.solver = self._get_solver()

    def _parse_input(self):
        """Convert input dict to class attrs."""
        for k, v in self.input_dict.items():
            # special handling for input path
            if k == "input_path":
                if not Path(v).exists():
                    raise InputFilesNotFoundError
                self.input_path = Path(v)
            else:
                setattr(self, k, v)

    def validate_input(self):
        """Validate against schema."""
        try:
            validate(self.input_dict, self.schema)
        except ValidationError as e:
            raise e
        return True

    def _get_solver(self):
        """Initialize solver object."""
        return solvers[self.solver_input["name"]](self.solver_input, self.input_path)

    def run(self):
        """Run convergence workflow."""
        # run at least once
        return_val = self.run_solver()
        self.results.append(return_val)
        while return_val < self.tol:
            return_val = self.run_solver()
            self.results.append(return_val)
        pass

    def run_solver(self) -> float:
        """Run solver and extract target converence value."""
        self.solver.run()
        return 0.0


class InputFilesNotFoundError(Exception):
    """Error when the provided input path string is not a valid Path."""
