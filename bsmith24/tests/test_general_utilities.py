import unittest
from unittest import mock

import sys
sys.path.append('..')
from src.general_utilities import General_Utilities


class Test_General_Utilities(unittest.TestCase):

    def setUp(self):
        self.template_file_path = '../template_files' 
 
    @mock.patch('builtins.open')
    def test_get_job_template(self, mock_builtin_open):
        General_Utilities.get_job_template(self.template_file_path) 
        builtin_open_calls = [mock.call(self.template_file_path+'/run.sh')]
        mock_builtin_open.assert_has_calls(builtin_open_calls)


if __name__ == '__main__':
   unittest.main()
