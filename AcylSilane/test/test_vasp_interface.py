import os
import shutil
import unittest

import ase.io
import numpy as np
from parameterized import parameterized

from AcylSilane.convtrack.calculation import Calculation
from AcylSilane.convtrack.convergence import Convergence


class local_io_tests(unittest.TestCase):
    def setUp(self):
        # Build Calculation object
        self.local_path = os.path.dirname(__file__)
        self.data_dir = os.path.join(os.path.dirname(self.local_path), "data/cu_vasp_files")
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.potcar = os.path.join(self.data_dir, "POTCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")

        # Build Convergence object
        self.output_path = os.path.join(self.local_path, "temp_io")
        # Ensure our temp directory doesn't exist already
        assert not os.path.isdir(self.output_path)
        self.convergence = Convergence(self.incar, self.poscar, self.potcar, max_size=3, root_dir=self.output_path,
                                       kpoints=self.kpoints, uniform_supercell=False)
        self.ase_obj = ase.io.read(self.poscar)

    def tearDown(self):
        if os.path.isdir(self.output_path):
            shutil.rmtree(self.output_path)

    def test_convergence_creates_parent_dir(self):
        self.assertTrue(os.path.isdir(self.output_path))

    @parameterized.expand([
        # Supercell size
        (2, 2, 2),
        (1, 1, 2),
        (1, 2, 1),
        (2, 1, 1),
        (3, 3, 3)
    ])
    def test_calculators_create_subfolders(self, a, b, c):
        subfolder = "".join(map(str, (a, b, c)))
        expected_path = os.path.join(self.output_path, subfolder)
        self.assertTrue(os.path.isdir(expected_path))

    def test_calculators_setup_calculations(self):
        cell_size = "222"
        expected_path = os.path.join(self.output_path, cell_size)
        self.convergence.calculations[cell_size].setup_calc()
        for filename in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            path_to_file = os.path.join(expected_path, filename)
            self.assertTrue(os.path.isfile(path_to_file))


class vasp_interface_tests(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(os.getcwd()), "data/sample_si_calc")
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.potcar = os.path.join(self.data_dir, "POTCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        self.calc = Calculation(self.incar, self.poscar, self.potcar,
                                dims=(1, 1, 1),
                                calc_folder=self.data_dir,
                                kpoints=self.kpoints)

    def tearDown(self):
        pass

    def test_calculation_contains_vasp_files(self):
        self.assertEqual(self.calc.incar, self.incar)
        self.assertEqual(self.calc.potcar, self.potcar)
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
        self.local_path = os.path.dirname(__file__)
        self.data_dir = os.path.join(os.path.dirname(self.local_path), "data/cu_vasp_files")
        self.output_path = os.path.join(self.local_path, "temp_io")
        # Ensure our temp directory doesn't exist already
        assert not os.path.isdir(self.output_path)
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.potcar = os.path.join(self.data_dir, "POTCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        self.calc = Calculation(self.incar, self.poscar, self.potcar,
                                dims=(1, 1, 1),
                                calc_folder=self.data_dir,
                                kpoints=self.kpoints)

        # Build Convergence object

        self.convergence = Convergence(self.incar, self.poscar, self.potcar, max_size=3, root_dir=self.output_path,
                                       kpoints=self.kpoints, uniform_supercell=False)
        self.ase_obj = ase.io.read(self.poscar)

    def tearDown(self):
        if os.path.isdir(self.output_path):
            shutil.rmtree(self.output_path)

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
    def test_convergence_creates_XxYxZ_calculation(self, a, b, c):
        dims = self.ase_obj.cell * (a, b, c)
        dict_key = "".join(map(str, (a, b, c)))
        self.assertTrue(np.all(self.convergence.calculations[dict_key].crystal.cell == dims))


class uniform_supercell_tests(unittest.TestCase):
    def setUp(self):
        # Build Calculation object
        self.local_path = os.path.dirname(__file__)
        self.data_dir = os.path.join(os.path.dirname(self.local_path), "data/cu_vasp_files")
        self.output_path = os.path.join(self.local_path, "temp_io")
        # Ensure our temp directory doesn't exist already
        assert not os.path.isdir(self.output_path)
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.potcar = os.path.join(self.data_dir, "POTCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        # Build Convergence object
        self.convergence = Convergence(self.incar, self.poscar, self.potcar, max_size=5, root_dir=self.output_path,
                                       kpoints=self.kpoints)
        self.ase_obj = ase.io.read(self.poscar)

    def tearDown(self):
        if os.path.isdir(self.output_path):
            shutil.rmtree(self.output_path)

    @parameterized.expand([
        # Supercell size
        (1, 1, 1),
        (2, 2, 2),
        (3, 3, 3),
        (4, 4, 4),
        (5, 5, 5)
    ])
    def test_convergence_uniform_creates_NxNxN_calculations(self, a, b, c):
        dims = self.ase_obj.cell * (a, b, c)
        dict_key = "".join(map(str, (a, b, c)))
        self.assertTrue(np.all(self.convergence.calculations[dict_key].crystal.cell == dims))

    @parameterized.expand([
        # Supercell size
        (1, 1, 2),
        (2, 2, 1),
        (3, 1, 1),
        (2, 4, 5),
        (5, 1, 3)
    ])
    def test_convergence_uniform_spacing_avoids_uneven_supercells(self, a, b, c):
        dict_key = "".join(map(str, (a, b, c)))
        self.assertNotIn(dict_key, self.convergence.calculations.keys())


if __name__ == "__main__":
    unittest.main()
