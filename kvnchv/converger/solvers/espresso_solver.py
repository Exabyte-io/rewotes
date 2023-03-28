"""Solver implementations for Quantum Espresso."""
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .base_solver import BaseSolver


class EspressoSolver(BaseSolver):
    """Derived class for quantum espresso simulations."""

    def __init__(self, input_dict: dict, input_path: Path, parameter_set: dict):
        super().__init__(input_dict, input_path, parameter_set)
        self.supported_parameters: list = ["k"]
        self.default_input_name: str = "pw.in"
        self._validate_parameters()
        self.params_string = "_".join(f"{k}-{v}" for k, v in self.parameter_set.items())

    def run(self):
        """Run solver subprocess."""
        print(f"can run on {self.solver_path}")

        # modify input file based on parameters
        run_input_file = f"{self.default_input_name}_{self.params_string}"
        run_output_dir = f"outdir_{self.params_string}"
        self._generate_new_file(run_input_file, run_output_dir)

        # run solver in input_path
        run_cmd = f"{self.solver_path} < {self.input_path.joinpath(run_input_file)}"
        result = subprocess.run(run_cmd, check=False, shell=True,
                                cwd=self.input_path, capture_output=True, text=True)
        output = result.stdout
        if result.returncode:
            raise SolverSubprocessFailedError(self._process_errors(output))

        self.parse_output(run_output_dir)

    def parse_output(self, outdir: Path):
        """Parse XML output and save as dict."""
        result_xml = self.input_path.joinpath(outdir, "__prefix__.xml")
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

    def _generate_new_file(self, new_input_fname: str, outdir: str = ''):
        """Modify input file to update specified parameters."""
        with open(self.input_path.joinpath(self.default_input_name), 'r') as default_fp:
            original_lines = default_fp.read()

        # set output dir
        new_lines = original_lines.replace("'outdir'", f"'{outdir}'")

        for param, value in self.parameter_set.items():
            regex, replace_string = getattr(self, f"_match_{param}")(value)
            new_lines = re.sub(regex, replace_string, new_lines)

        with open(self.input_path.joinpath(new_input_fname), 'w+') as new_fp:
            new_fp.write(new_lines)

    def _match_k(self, k: int):
        """Return match regex and replace_string.

        Room for flexibility improvement. Assumptions:
        - K_POINTS automatic is used.
        - K_POINTS is the last command issued
        """
        return ("K_POINTS automatic.*",
                f"K_POINTS automatic\n{k} {k} {k} 0 0 0")


class SolverSubprocessFailedError(Exception):
    """Return processed subprocess error."""
