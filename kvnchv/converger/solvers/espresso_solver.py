"""Solver implementations for Quantum Espresso."""
import re
import subprocess
from pathlib import Path

from .base_solver import BaseSolver


class EspressoSolver(BaseSolver):
    """Derived class for quantum espresso simulations."""

    def __init__(self, input_dict: dict, input_path: Path, parameter_set):
        super().__init__(input_dict, input_path, parameter_set)
        self.supported_parameters = ["k"]
        self._validate_parameters()

    def run(self):
        """Excute solver subprocess."""
        print(f"can run on {self.solver_path}")
        # run solver in input_path
        run_cmd = f"{self.solver_path} < {self.input_path.joinpath('pw.in')}"
        result = subprocess.run(run_cmd, check=False, shell=True,
                                cwd=self.input_path, capture_output=True, text=True)
        output = result.stdout
        if result.returncode:
            raise SolverSubprocessFailedError(self._process_errors(output))

    def parse_output(self):
        """Parse XML output and save as dict."""
        """
        result = {
            "val1": 123,
            "val2": 456,
        } #validate against OneOf type schema?
        """
        pass

    @staticmethod
    def _process_errors(stdout: str):
        """Identify exact error for failed subprocess."""
        # espresso does not output CRASH in the cwd
        return re.search(r"Error.*\s*.+", stdout)[0]

    def _search_and_replace_params(self):
        pass


class SolverSubprocessFailedError(Exception):
    """Return processed subprocess error."""
