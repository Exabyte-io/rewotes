import itertools
import os
import typing

import ase.build
import ase.io


class Calculation(object):
    def __init__(self, incar: str, poscar: str, potcar: str,
                 dims: typing.Iterable, calc_folder: str, kpoints: str = None):
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
        self.potcar = potcar
        self.kpoints = kpoints

        self.calc_folder = calc_folder
        if not os.path.isdir(self.calc_folder):
            os.mkdir(self.calc_folder)

    def setup_calc(self):
        to_copy = {"INCAR" : self.incar,
                   "KPOINTS" : self.kpoints,
                   "POTCAR" : self.potcar}
        for filename, source in to_copy.items():
            path_to_file = os.path.join(self.calc_folder, filename)
            if source is not None:
                assert not os.path.isfile(path_to_file)
                with open(source, "r") as inp, open(path_to_file, "w") as outp:
                    for line in inp:
                        outp.write(line)
        ase.io.write(os.path.join(self.calc_folder, "POSCAR"), self.crystal, format="vasp")


    def started(self):
        """
        Checks whether a calculation has started or not
        :return: True if the job has started, otherwise False
        """
        contents = os.listdir(self.calc_folder)
        if "OUTCAR" in contents:
            started = True
        else:
            started = False
        return started

    def complete(self):
        """
        Checks whether the corresponding VASP job has completed
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

    def energy(self):
        """
        Checks the OUTCAR for energy
        :return: A float containing the energy if the job has completed, otherwise None
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



class Convergence(object):
    def __init__(self, incar, poscar, potcar, max_size, root_dir, kpoints=None):
        """
        Automatically tracks convergence and generates new calculations.

        :param incar: Path to the INCAR
        :type incar: str
        :param poscar: Path to the POSCAR
        :type poscar: str
        :param potcar: Path to the POTCAR
        :type potcar: str
        :param max_size: Maximum supercell size to consider
        :type max_size: int
        :param root_dir: Root directory for the convergence calculation, locally
        :type root_dir: str
        :param kpoints: Path to the KPOINTS
        :type kpoints: str
        """
        self.incar = incar
        self.poscar = poscar
        self.potcar = potcar
        self.kpoints = kpoints
        self.max_size = max_size
        self.created_dir = False
        if root_dir is not None:
            self.create_directory_skeleton(root_dir)

        self.calculations = {}

        # Create the calculations and assign subfolders
        cell_sizes = list(itertools.product(range(1, max_size + 1),
                                            range(1, max_size + 1),
                                            range(1, max_size + 1)
                                            )
                          )
        for size in cell_sizes:
            dims = "".join(map(str, size))
            foldname = os.path.join(root_dir, dims)
            calculation = Calculation(self.incar, self.poscar,self.potcar, size, foldname, self.kpoints)
            self.calculations[dims] = calculation

    def create_directory_skeleton(self, root_dir):
        # Make sure we haven't made the directory structure before
        assert not self.created_dir
        self.root_dir = root_dir
        self.created_dir = True
        # Make the parent directory
        os.mkdir(root_dir)



if __name__ == "__main__":
    pass
