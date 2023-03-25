"""Solver implementations for Quantum Espresso."""
import subprocess
from pathlib import Path

from .base_solver import BaseSolver


class EspressoSolver(BaseSolver):
    """Derived class for quantum espresso simulations."""

    def __init__(self, input_dict: dict, input_path: Path):
        super().__init__(input_dict, input_path)

    def run(self):
        """Excute solver subprocess."""
        print(f"can run on {self.solver_path}")
        pass

    def parse_output(self):
        """Parse XML output and save as dict."""
        """
        result = {
            "val1": 123,
            "val2": 456,
        } #validate against OneOf type schema?
        """
        pass
