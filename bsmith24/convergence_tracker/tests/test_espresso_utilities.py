import unittest
from unittest import mock

import sys
sys.path.append('..')
from espresso_utilities import Espresso_Calculation

class Test_Espresso_Utilities(unittest.TestCase):

    def setUp(self):
        self.template_file_path = '../template_files' 
        self.ecutwfcs = ['26','28','30', '32']
        self.kpoints = ['4 4 4 0 0 0']
        self.espresso_calculation = Espresso_Calculation(self.template_file_path, self.ecutwfcs, self.kpoints)
        self.fake_readlines = ['foo\n', 'bar\n', 'foobar\n']

    def test_espresso_calculation_initialization(self):
        self.assertEqual(self.espresso_calculation.template_file_path, self.template_file_path)
        self.assertEqual(self.espresso_calculation.ecutwfcs, self.ecutwfcs)
        self.assertEqual(self.espresso_calculation.kpoints, self.kpoints)

    @mock.patch('builtins.open', new_callable=mock.mock_open, read_data='foo\nbar\nfoobar\n')
    def test_get_espresso_job_template(self, mock_builtin_open):
        self.espresso_calculation.get_espresso_job_template()
        builtin_open_calls = [mock.call(self.template_file_path+'/run.sh')]
        mock_builtin_open.assert_has_calls(builtin_open_calls)
        self.assertEqual(self.espresso_calculation.espresso_job_template, self.fake_readlines)

    def test_update_espresso_job_template_AssertionError(self):
        self.assertRaises(AssertionError, self.espresso_calculation.update_espresso_job_template)

    def test_update_espresso_job_template_assignments(self):
        self.espresso_calculation.espresso_job_template = self.fake_readlines
        self.espresso_calculation.update_espresso_job_template()
        self.assertEqual(self.espresso_calculation.espresso_job_template[1], 'kpoints='+"'"+' '.join(self.kpoints)+"'"+'\n') 
        self.assertEqual(self.espresso_calculation.espresso_job_template[2], 'for ecut in '+' '.join(self.ecutwfcs)+'\n')

    @mock.patch('subprocess.getoutput')
    def test_get_espresso_total_energy(self, mock_subprocess_getoutput):
        mock_subprocess_getoutput.return_value = '! fake total energy 1.0 Ry'
        fake_total_energy = self.espresso_calculation.get_espresso_total_energy('fake_espresso_output.out')
        assert fake_total_energy == 13605.662285137  

    @mock.patch('subprocess.getoutput')
    def test_get_each_espresso_total_energy(self, mock_subprocess_getoutput):
        meV2Ry = 1.0/13605.662285137
        mock_iterables = ['! fake total energy '+str(i*meV2Ry)+' Ry' for i in range(len(self.ecutwfcs))]
        mock_subprocess_getoutput.side_effect = mock_iterables
        self.espresso_calculation.get_each_espresso_total_energy()
        assert self.espresso_calculation.total_energies == [0.0, 1.0, 2.0, 3.0]


if __name__ == '__main__':
   unittest.main()
