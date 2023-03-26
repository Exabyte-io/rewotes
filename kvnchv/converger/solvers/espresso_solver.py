"""Solver implementations for Quantum Espresso."""
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

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

        self.parse_output()

    def parse_output(self):
        """Parse XML output and save as dict."""
        result_xml = self.input_path.joinpath("outdir", "__prefix__.xml")
        tree = ET.parse(result_xml)
        self.results["etot"] = self.parse_etot(tree)
        self.results["fnorm"] = self.parse_forces(tree)

    @staticmethod
    def parse_etot(tree: ET) -> float:
        """Get value of etot in result tree."""
        return float(tree.findall('output/total_energy/etot')[0].text)

    @staticmethod
    def parse_forces(tree: ET) -> float:
        """Return norm of force vector on all atoms."""
        forces = tree.findall('output/forces')[0].text.strip().split("\n")
        fvecs = np.array([f.split() for f in forces], dtype='float')
        return np.linalg.norm(fvecs)

    @staticmethod
    def _process_errors(stdout: str):
        """Identify exact error for failed subprocess."""
        # espresso does not output CRASH in the cwd
        return re.search(r"Error.*\s*.+", stdout)[0]

    def _search_and_replace_params(self):
        pass


class SolverSubprocessFailedError(Exception):
    """Return processed subprocess error."""
