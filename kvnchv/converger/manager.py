"""Convergence workflow manager."""
import itertools
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
from jsonschema import ValidationError, validate

from .schema import input_schema
from .solvers import EspressoSolver

solvers = {"espresso": EspressoSolver}


class Manager():
    """Main class for convergence workflow run."""

    def __init__(self, input_dict: dict):
        self.input_dict: dict = deepcopy(input_dict)
        self.schema: dict = input_schema  # schema to validate against

        self.solver_input: dict  # solver-configuration
        self.input_path: Path
        self.target: str  # convergence target value name
        self.tol: float  # tolerance (absolute) to satisfy
        self.parameter_space: list  # list of input parameter specifications

        self.data: pd.DataFrame  # DataFrame for parameter sets and outputs

        self.validate_input()
        self._parse_input()

    def _parse_input(self):
        """Convert input dict to class attrs."""
        for param, val in self.input_dict.items():
            # special handling for input path
            if param == "input_path":
                if not Path(val).exists():
                    raise InputFilesNotFoundError
                self.input_path = Path(val)
            else:
                setattr(self, param, val)

    def validate_input(self):
        """Validate against schema."""
        try:
            validate(self.input_dict, self.schema)
        except ValidationError as err:
            raise err
        return True

    def _get_solver(self, run_params: dict):
        """Initialize solver object."""
        return solvers[self.solver_input["name"]](
                    self.solver_input,
                    self.input_path,
                    run_params
                )

    def run(self):
        """Run convergence workflow."""
        # run at least once

        self._process_parameters()

        target_eval = self.run_solver(0)
        run_idx = 1
        while target_eval < self.tol:
            target_eval = self.run_solver(run_idx)
            run_idx += 1

    def run_solver(self, idx: int) -> float:
        """Run solver for parameter set 'idx' and extract target converence value."""
        run_params = self.data.iloc[idx]
        solver = self._get_solver(run_params)
        solver.run()

        return 0.0

    def _process_parameters(self) -> pd.DataFrame:
        """Process parameter_space input values define simulation set."""
        param_arrays = []
        param_names = []
        for param_d in deepcopy(self.parameter_space):
            # if delta is not provided, assume unit step
            del_p = 1.0 if "delta" not in param_d else param_d["delta"]
            param_names.append(param_d["name"])
            param_arrays.append(np.arange(param_d["start"], param_d["max"], del_p))

        self.data = pd.DataFrame(
            list(itertools.product(*param_arrays)),
            columns=param_names)


class InputFilesNotFoundError(Exception):
    """Error when the provided input path string is not a valid Path."""
