import itertools
import os

from AcylSilane.convtrack.calculation import Calculation


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