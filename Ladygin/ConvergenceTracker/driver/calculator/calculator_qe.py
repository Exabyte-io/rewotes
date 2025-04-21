import subprocess as sp
import os

from .calculator_base import calculator
from .utils import wait_qe

class calculator_qe(calculator):
    """Class for point calculations using quantum espresso"""


    def __init__(self, path_to_exec:str = "/pscratch/sd/v/vladygin/tools/q-e_new/q-e/bin") -> None:
        """path_to_exec - path to executable"""
        self.path_to_exec = path_to_exec


    def run_scf(self, input_file:str, output_file:str) -> None:
        """Runs scf calculations
        
           input_file - input settings
           ouput_file - output file to save results
        """
        binary = os.path.join(self.path_to_exec, "pw.x")
        sp.run(f"{binary} -input {input_file} > {output_file}", shell = True)

class calculator_qe_par(calculator_qe):
    """Class for point calculations using quantum espresso using parralelization along q"""


    def __init__(self, path_to_exec:str = "/pscratch/sd/v/vladygin/tools/q-e_new/q-e/bin", 
                 ncores: int = 1, nk: int = 1,
                 job: bool = False,
                 workdir: str = None) -> None:
        """path_to_exec - path to executable (kept for testing perposes better remove)
           ncores - number of cores for parallelization
           nk - quantum espresso parallelization setting along q   
        """
        self.path_to_exec = path_to_exec

        self.ncores = ncores
        self.nk = nk

        self.job = job
        if job:
            if type(workdir) == type(None):
                raise ValueError("For job calculation, one should provide a directory")
            self.workdir = workdir

    def run_scf(self, input_file:str, output_file:str) -> None:
        """Runs scf calculations
        
           input_file - input settings
           ouput_file - output file to save results
           job - submit job or run in shell
        """
        binary = os.path.join(self.path_to_exec, "pw.x")

        if not self.job:
            sp.run(f"srun -n {self.ncores} {binary} -nk {self.nk} -input {input_file} > {output_file}", shell = True)
        else:
            job_file = os.path.join(self.workdir, "calc.sh")
            # Designed to run on one node
            if self.ncores > 128:
                raise ValueError("Number of cores is greater than number of cores per node")
                
            with open(f"{job_file}", "a") as file:
                file.write(f"#!/bin/bash\n")
                file.write(f"srun -n {self.ncores} {binary} -nk {self.nk} -input {input_file} > {output_file}\n")

            os.remove(output_file)
            sp.run(f"sbatch -J qe_scf -q premium -C cpu -t 1:00:00 -N 1 -n {self.ncores} {job_file}", shell = True)
            wait_qe(output_file)