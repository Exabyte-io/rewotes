import itertools
import os

from AcylSilane.convtrack.calculation import Calculation


class Convergence(object):
    def __init__(self, incar, poscar, potcar, max_size, root_dir, kpoints=None, uniform_supercell=True):
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
        :param uniform_supercell: Whether to only consider NxNxN supercells, or to allow uneven sizes
        :type uniform_supercell: bool
        """
        self.incar = incar
        self.poscar = poscar
        self.potcar = potcar
        self.kpoints = kpoints
        self.max_size = max_size
        self.uniform_supercell = uniform_supercell
        self.created_dir = False
        if root_dir is not None:
            self.create_directory_skeleton(root_dir)

        self.calculations = {}

        self.create_calculations()

    def create_directory_skeleton(self, root_dir):
        """
        Sets up a calculation in a root dir, setting self.root_dir if it hasn't been set before.
        :param root_dir:
        :return:
        """
        # Make sure we haven't made the directory structure before
        assert not self.created_dir
        self.root_dir = root_dir
        self.created_dir = True
        # Make the parent directory
        os.mkdir(root_dir)

    def create_calculations(self):
        """
        Creates the actual Calculation objects, and lets them create their subfolders.
        :return:
        """
        if self.uniform_supercell:
            # Create supercells with uniform repetition in each direction up to self.max_size
            cell_sizes = [(i,i,i) for i in range(1, self.max_size + 1)]
        else:
            # Create all possible supercells (e.g. [1,1,1], [1,1,2], [1,2,1], etc) up to self.max_size
            cell_sizes = list(itertools.product(range(1, self.max_size + 1),
                                                range(1, self.max_size + 1),
                                                range(1, self.max_size + 1)
                                                )
                              )
        for size in cell_sizes:
            dims = "".join(map(str, size))
            foldname = os.path.join(self.root_dir, dims)
            calculation = Calculation(self.incar, self.poscar, self.potcar, size, foldname, self.kpoints)
            self.calculations[dims] = calculation


if __name__ == "__main__":
    pass
