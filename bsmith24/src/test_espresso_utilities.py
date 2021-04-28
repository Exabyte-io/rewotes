import unittest
from unittest import mock

from espresso_utilities import Espresso_Calculation


class Test_Espresso_Utilities(unittest.TestCase):

    def setUp(self):
        self.template_files_path = '../files' 
        self.ecutwfcs = ['26','28','30', '32']
        self.kpoints = ['4 4 4 0 0 0']
        self.espresso_calculation = Espresso_Calculation(self.template_files_path, self.ecutwfcs, self.kpoints)

    def test_espresso_calculation_initialization(self):
        self.assertEqual(self.espresso_calculation.template_files_path, self.template_files_path)
        self.assertEqual(self.espresso_calculation.ecutwfcs, self.ecutwfcs)
        self.assertEqual(self.espresso_calculation.kpoints, self.kpoints)
 
    @mock.patch('builtins.open')
    def test_getting_espresso_input_template(self, mock_builtin_open):
        self.espresso_calculation.get_espresso_input_template() 
        builtin_open_calls = [mock.call('../files/pw.in')]
        mock_builtin_open.assert_has_calls(builtin_open_calls)

    @mock.patch('builtins.open')
    def test_getting_espresso_submit_template(self, mock_builtin_open):
        self.espresso_calculation.get_espresso_submit_template()
        builtin_open_calls = [mock.call('../files/job.pbs')]
        mock_builtin_open.assert_has_calls(builtin_open_calls)

if __name__ == '__main__':
   unittest.main()
