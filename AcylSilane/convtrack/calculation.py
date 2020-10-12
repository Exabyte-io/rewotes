import itertools
import typing

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
