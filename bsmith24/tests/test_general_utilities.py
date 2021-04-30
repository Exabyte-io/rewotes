import unittest
from unittest import mock

import sys
sys.path.append('..')
from src.general_utilities import General_Utilities


class Test_General_Utilities(unittest.TestCase):

    def setUp(self):
        self.template_file_path = '../template_files'
        self.fake_readlines = ['foo\n', 'bar\n', 'foobar\n'] 
 
    @mock.patch('builtins.open', new_callable=mock.mock_open, read_data='foo\nbar\nfoobar\n')
    def test_get_job_template(self, mock_builtin_open):
        filedata = General_Utilities.get_job_template(self.template_file_path) 
        builtin_open_calls = [mock.call(self.template_file_path+'/run.sh')]
        mock_builtin_open.assert_has_calls(builtin_open_calls)
        self.assertEqual(self.fake_readlines, filedata)

    @mock.patch('builtins.open')
    def test_write_driver_script(self, mock_builtin_open):
        fake_list_of_strings = ['foo\n', 'bar\n', 'foobar\n']
        General_Utilities.write_driver_script(fake_list_of_strings)
        builtin_open_calls = [mock.call(string_var) for string_var in fake_list_of_strings]
        mock_builtin_open.return_value.__enter__().write.assert_has_calls(builtin_open_calls)

    # Just setting up a few test cases for illustration, but could use @pytest.mark.parameterize ...
    def test_is_converged_1(self):
        values, tolerance = [10.0,7.0,5.0,4.0,3.5,3.25], 1.0
        self.assertEqual(General_Utilities.is_converged(values, tolerance), [True, 5.0, 2])

    def test_is_converged_2(self):
        values, tolerance = [10.0,7.0,5.0,4.0,3.5,3.25], 0.5
        self.assertEqual(General_Utilities.is_converged(values, tolerance), [True, 4.0, 3])

    def test_is_converged_3(self):
        values, tolerance = [10.0,7.0,5.0,4.0,3.5,3.25], 0.1
        self.assertEqual(General_Utilities.is_converged(values, tolerance), [False, None, None])

if __name__ == '__main__':
   unittest.main()
