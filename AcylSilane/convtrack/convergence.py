import itertools
import os
import typing

from .calculation import Calculation


class Convergence(object):
    def __init__(self, incar: str, poscar: str, potcar: str, max_size: int, root_dir: str,
                 kpoints: str = None, uniform_supercell: bool = True):
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


    def create_directory_skeleton(self, root_dir: str) -> None:
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

    def create_calculations(self) -> None:
        """
        Creates the actual Calculation objects, and lets them create their subfolders.
        :return:
        """
        if self.uniform_supercell:
            # Create supercells with uniform repetition in each direction up to self.max_size
            cell_sizes = [(i, i, i) for i in range(1, self.max_size + 1)]
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
            calculation = Calculation(self.incar, self.poscar, self.potcar, size, foldname, self.kpoints,
                                      create=True)
            self.calculations[dims] = calculation

    def submit_calculations(self) -> None:
        """
        Submits the calculations
        :return: None
        """
        for key, calculation in self.calculations.items():
            calculation.submit()

    def all_calculations_done(self) -> bool:
        """
        Returns the status of all calculations being held in the Convergence object.
        :return: True if all are done, else false.
        """
        for cellsize, calculation in self.calculations.items():
            if not calculation.complete():
                return False
        return True

    def energies_to_csv(self, filename: str) -> None:
        """
        Writes the calculation energies to a CSV
        """
        energies = {}
        if self.all_calculations_done():
            for cellsize, calculation in self.calculations.items():
                energies[cellsize] = calculation.energy()
        with open(filename, "w") as outp:
            outp.write("Cell_Dimensions,Energy_eV\n")
            for cellsize, energy in energies:
                outp.write(f"{cellsize},{energy}\n")

    def calc_converged_size(self, convergence_criterion:float) -> typing.List[int]:
        """
        Returns the supercell size at which the energy has converged
        :param convergence_criterion: Convergence criterion, in eV/atom
        :return:
        """
        if len(self.calculations) <= 1:
            raise ValueError("Too few calculations to determine convergence. You need at least two.")

        # Todo: Extend this method to non-uniform cells
        if not self.uniform_supercell:
            raise NotImplementedError("Current version only supports uniformly-sized supercells.")

        ref_energy = self.calculations["111"].intensive_energy()
        best_convergence = None
        best_i = None
        for i in range(2,self.max_size+1):
            # We're checking the energy difference relative to the next-smallest calculation
            cell_key = f"{i}{i}{i}"
            abs_energy = self.calculations[cell_key].intensive_energy()
            delta_energy = abs(abs_energy - ref_energy)

            # Log best convergence, in case we don't actually satisfy the convergence criterion
            if best_convergence is None:
                best_convergence = delta_energy
                best_i = i
            elif delta_energy < best_convergence:
                best_convergence = delta_energy
                best_i = i

            # Break early if we're converged in supercell size
            if delta_energy < convergence_criterion:
                return [i,i,i]
            ref_energy = abs_energy

        raise ValueError(f"Supercells are not converged. Best convergence is {best_convergence} eV at {best_i}/{best_i}/{best_i}.")







if __name__ == "__main__":
    pass
