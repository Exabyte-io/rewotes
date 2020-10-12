import os
import unittest

import ase.io
import numpy as np
from parameterized import parameterized

from AcylSilane.convtrack.calculation import Calculation, Convergence


class vasp_interface_tests(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(os.getcwd()), "data/sample_si_calc")
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        self.calc = Calculation(self.incar, self.poscar,
                                dims=(1, 1, 1),
                                calc_folder=self.data_dir,
                                kpoints=self.kpoints)

    def tearDown(self):
        pass

    def test_calculation_contains_vasp_files(self):
        self.assertEqual(self.calc.incar, self.incar)
        self.assertEqual(self.calc.poscar, self.poscar)
        self.assertEqual(self.calc.kpoints, self.kpoints)

    def test_detects_job_has_started(self):
        self.assertTrue(self.calc.started())

    def test_detects_job_has_not_started(self):
        self.calc.calc_folder = "../data/cu_vasp_files"
        self.assertFalse(self.calc.started())

    def test_calculation_detects_finished_job(self):
        self.assertTrue(self.calc.complete())

    def test_calculation_returns_energy_if_complete(self):
        true_energy = -8.31258263
        self.assertEqual(self.calc.energy(), true_energy)

    def test_calculation_returns_None_energy_if_incomplete(self):
        self.calc.calc_folder = "../data/cu_vasp_files"
        self.assertIsNone(self.calc.energy())


class supercell_creation_tests(unittest.TestCase):
    def setUp(self):
        # Build Calculation object
        self.data_dir = os.path.join(os.path.dirname(os.getcwd()), "data/cu_vasp_files")
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        self.calc = Calculation(self.incar, self.poscar,
                                dims=(1, 1, 1),
                                calc_folder=self.data_dir,
                                kpoints=self.kpoints)

        # Build Convergence object
        self.convergence = Convergence(self.incar, self.poscar, max_size=3, kpoints=self.kpoints)
        self.ase_obj = ase.io.read(self.poscar)

    def tearDown(self):
        pass

    def test_convergence_creates_1x1x1_calculation(self):
        self.assertTrue(np.all(self.convergence.calculations["111"].crystal.cell == self.ase_obj.cell))

    @parameterized.expand([
        # Supercell size
        (2, 2, 2),
        (1, 1, 2),
        (1, 2, 1),
        (2, 1, 1),
        (3, 3, 3)
    ])
    def test_convergence_creates_NxNxN_calculation(self, a, b, c):
        dims = self.ase_obj.cell * (a, b, c)
        dict_key = "".join(map(str, (a, b, c)))
        self.assertTrue(np.all(self.convergence.calculations[dict_key].crystal.cell == dims))


if __name__ == "__main__":
    unittest.main()
