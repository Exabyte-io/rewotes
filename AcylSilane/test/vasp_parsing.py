import unittest
import os

from AcylSilane.convtrack.calculation import Calculation


class CalculationTests(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(os.getcwd()), "data")
        self.incar = os.path.join(self.data_dir, "INCAR")
        self.poscar = os.path.join(self.data_dir, "POSCAR")
        self.kpoints = os.path.join(self.data_dir, "KPOINTS")
        self.calc = Calculation(self.incar, self.poscar, kpoints=self.kpoints)

    def tearDown(self):
        pass

    def test_calculation_contains_vasp_files(self):
        self.assertEqual(self.calc.incar, self.incar)
        self.assertEqual(self.calc.poscar, self.poscar)
        self.assertEqual(self.calc.kpoints, self.kpoints)


if __name__ == "__main__":
    unittest.main()