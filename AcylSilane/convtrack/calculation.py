import itertools
import typing
import os
import ase.io, ase.build


class Calculation(object):
    def __init__(self, incar: str, poscar: str, dims: typing.Iterable, kpoints: str = None):
        """
        Object containing methods related to a VASP calculation.

        :param incar: Path to the INCAR
        :type incar: str
        :param poscar: Path to the POSCAR
        :type poscar: str
        :param dims: An iterable containing the unit cell's A/B/C size, in terms of number of repeating units.
        :type dims: typing.Iterable
        :param kpoints: Path to the KPOINTS
        :type kpoints: str
        """
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints
        self.crystal = ase.io.read(self.poscar)*dims

        self.calc_folder = "../data/sample_si_calc"

    def started(self):
        contents = os.listdir(self.calc_folder)
        if "OUTCAR" in contents:
            started = True
        else:
            started = False
        return started

    def complete(self):
        """
        Checks whether the corresponding VASP job has completed
        :return:
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



class Convergence(object):
    def __init__(self, incar, poscar, max_size, kpoints=None):
        """
        Automatically tracks convergence and generates new calculations.

        :param incar: Path to the INCAR
        :type incar: str
        :param poscar: Path to the POSCAR
        :type poscar: str
        :param kpoints: Path to the KPOINTS
        :type kpoints: str
        """
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints
        self.max_size = max_size

        self.calculations = {}
        cell_sizes = list(itertools.product(range(1, max_size + 1),
                                            range(1, max_size + 1),
                                            range(1, max_size + 1)
                                            )
                          )
        for size in cell_sizes:
            dims = "".join(map(str,size))
            calculation = Calculation(self.incar, self.poscar, size, self.kpoints)
            self.calculations[dims] = calculation


if __name__ == "__main__":
    pass
