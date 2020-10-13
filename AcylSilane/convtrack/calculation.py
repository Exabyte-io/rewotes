import os
import typing

import ase.build
import ase.io

from .queueinterface import QueueInterface


class Calculation(object):
    def __init__(self, incar: str, poscar: str, potcar: str,
                 dims: typing.Iterable, calc_folder: str, kpoints: str = None,
                 create=False, submit=False):
        """
        Object containing methods related to a VASP calculation.

        :param incar: Path to the INCAR
        :type incar: str
        :param poscar: Path to the POSCAR
        :type poscar: str
        :param potcar: Path to the POTCAR
        :type potcar: str
        :param dims: An iterable containing the unit cell's A/B/C size, in terms of number of repeating units.
        :type dims: typing.Iterable
        :param calc_folder: Location on the drive for the calculation
        :type calc_folder: str
        :param kpoints: Path to the KPOINTS
        :type kpoints: str
        """
        self.incar = incar
        self.crystal = ase.io.read(poscar) * dims
        self.dims = dims
        self.potcar = potcar
        self.kpoints = kpoints
        self.job_id = None
        self.calc_folder = os.path.abspath(calc_folder)
        if not os.path.isdir(self.calc_folder):
            os.mkdir(self.calc_folder)
        if create:
            self.setup_calc(submit)

    def setup_calc(self, submit=False) -> None:
        """
        Sets up a calculation, by copying in the INCAR / KPOINTS / POTCAR, and writing a POSCAR.
        :return:  None
        """
        to_copy = {"INCAR": self.incar,
                   "KPOINTS": self.kpoints,
                   "POTCAR": self.potcar}
        for filename, source in to_copy.items():
            path_to_file = os.path.join(self.calc_folder, filename)
            if source is not None:
                assert not os.path.isfile(path_to_file)
                with open(source, "r") as inp, open(path_to_file, "w") as outp:
                    for line in inp:
                        outp.write(line)
        ase.io.write(os.path.join(self.calc_folder, "POSCAR"), self.crystal, format="vasp")

        # Todo: make calc parameters non non-hardcoded
        script_location = os.path.join(self.calc_folder, "job_vasp.pbs")
        QueueInterface.write_submission_script(script_location,
                                               jobname=f"Cu_{''.join(list(map(str, self.dims)))}",
                                               n_nodes=1,
                                               cores_per_node=8,
                                               walltime="00:00:30:00",
                                               quality_of_service="OR",
                                               email="jrd101@pitt.edu")
        if submit:
            cwd = os.getcwd()
            os.chdir(self.calc_folder)
            jobid = QueueInterface.qsub("job_vasp.pbs")
            os.chdir(cwd)

    def started(self) -> bool:
        """
        Checks whether a calculation has started or not by looking at whether OUTCAR has been created.
        :return: True if the job has started, otherwise False
        """
        contents = os.listdir(self.calc_folder)
        if "OUTCAR" in contents:
            started = True
        else:
            started = False
        return started

    def complete(self) -> bool:
        """
        Checks whether the corresponding VASP job has completed by checking whether the timing info is in the OUTCAR.
        :return: True if the job has completed, otherwise False
        """
        complete = False
        if self.started():
            # The timing information is only reported when the job successfully ends, and as far as I know is the only
            # consistent/reliable method to determine whether the code ended without error.
            termination_indicator = "General timing and accounting informations for this job"
            with open(os.path.join(self.calc_folder, "OUTCAR"), "r") as outcar:
                for line in outcar:
                    if termination_indicator in line:
                        complete = True
                        break
        return complete

    def energy(self) -> float:
        """
        Checks the OUTCAR for energy
        :return: A float containing the energy (in eV) if the job has completed, otherwise None
        """
        energy = None
        if self.complete():
            # Specifically, we want this one because some smearing methods need to be extrapolated to sigma=0
            indicator = "energy(sigma->0)"
            with open(os.path.join(self.calc_folder, "OUTCAR"), "r") as outcar:
                for line in outcar:
                    if indicator in line:
                        # Energy is the last value on the line
                        energy = float(line.strip().split()[-1])
        return energy


if __name__ == "__main__":
    pass
