import unittest
from unittest import mock

import sys
sys.path.append('..')
from src.espresso_utilities import Espresso_Calculation


class Test_Espresso_Utilities(unittest.TestCase):

    def setUp(self):
        self.template_file_path = '../template_files' 
        self.ecutwfcs = ['26','28','30', '32']
        self.kpoints = ['4 4 4 0 0 0']
        self.espresso_calculation = Espresso_Calculation(self.template_file_path, self.ecutwfcs, self.kpoints)

    def test_espresso_calculation_initialization(self):
        self.assertEqual(self.espresso_calculation.template_file_path, self.template_file_path)
        self.assertEqual(self.espresso_calculation.ecutwfcs, self.ecutwfcs)
        self.assertEqual(self.espresso_calculation.kpoints, self.kpoints)

    def test_update_espresso_job_template_assert(self):
        self.assertRaises(AssertionError, self.espresso_calculation.update_espresso_job_template)

    def test_update_espresso_job_template_assignments(self):
        fake_list_of_strings = ['foo\n', 'bar\n', 'foobar\n']
        self.espresso_calculation.espresso_job_template = fake_list_of_strings
        self.espresso_calculation.update_espresso_job_template()
        self.assertEqual(self.espresso_calculation.espresso_job_template[1], 'kpoints='+"'"+' '.join(self.kpoints)+"'"+'\n') 
        self.assertEqual(self.espresso_calculation.espresso_job_template[2], 'for ecut in '+' '.join(self.ecutwfcs)+'\n')


if __name__ == '__main__':
   unittest.main()
