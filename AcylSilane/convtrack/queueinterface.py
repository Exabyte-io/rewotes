import os
import re
import subprocess
import typing


class QueueInterface(object):
    def __init__(self):
        pass

    @staticmethod
    def write_submission_script(path: str,
                                jobname: str, n_nodes: int, cores_per_node: int, walltime: str,
                                quality_of_service: str, email: str = None) -> None:
        """
        Writes a submission script to a specified path. Will raise a FileExistsError if a file already exists there.
        :param path: Path to the submission script.
        :type path: str
        :param jobname: Name of the job.
        :type jobname: str
        :param n_nodes:  Number of nodes.
        :type n_nodes: int
        :param cores_per_node: Number of cores per node.
        :type cores_per_node: int
        :param walltime: Walltime, as a string, in a format acceptable for the PBS scheduler.
        :type walltime: str
        :param quality_of_service: Quality of service string, per the Exabyte documentation. Examples include "D", "OR", "OR4", etc.
        :type quality_of_service: str
        :param email: Optional. Which email to contact when a job completes.
        :type email: str
        :return: None
        """
        if os.path.isfile(path):
            raise FileExistsError(f"File at {path} already exists.")
        lines = ["#!/bin/bash",
                 f"#PBS -N {jobname}",
                 "#PBS -j oe",
                 f"#PBS -l nodes={n_nodes}",
                 f"#PBS -l ppn={cores_per_node}",
                 f"#PBS -l walltime={walltime}",
                 f"#PBS -q {quality_of_service}"]
        if email:
            lines += ["#PBS -m abe",
                      f"#PBS -M {email}"]
        lines += ["module add vasp/535-i-174-impi-044",
                  "cd $PBS_O_WORKDIR",
                  "mpirun -np $PBS_NP vasp > vasp.log"]
        with open(path, "w") as script:
            script.write("\n".join(lines))

    @staticmethod
    def qstat() -> typing.List[typing.Dict]:
        """
        Calls the system's qstat and parses its results.
        :return: A list of dictionaries, each containing job information.
        """
        result = subprocess.check_output("qstat")
        if isinstance(result, bytes):
            result = result.decode("utf-8")
        raw_qstat = result.split("\n")

        # Format is Headers, followed by a series of dashed lines, followed by actual data
        keys = re.split("\s{2,}", raw_qstat[0])
        vals = []
        for index, line in enumerate(raw_qstat):
            if index <= 1:
                # Skip header and decorative dashes
                continue
            elif re.match("^\s*$", line):
                # Skip blank lines
                continue
            else:
                vals.append(re.split("\s{2,}", line.strip()))
        return [dict(zip(keys, val)) for val in vals]

    @staticmethod
    def qsub(path_to_file: str) -> str:
        """
        Calls the system's qsub method, and returns a job ID for tracking purposes.
        :param path_to_file: path to the file being called
        :type path_to_file: str
        :return: None
        """
        job_id = subprocess.check_output(['qsub', path_to_file])
        if isinstance(job_id, bytes):
            job_id.decode("utf-8")
        return job_id
