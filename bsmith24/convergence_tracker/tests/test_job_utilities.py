import unittest

import sys
sys.path.append('..')
from job_utilities import Job
from job_utilities import Submit_Utilities


class Test_Submit_Utilities(unittest.TestCase):

    def test_submit_initialization(self):
        submit_tool = Submit_Utilities()

class Test_Job(unittest.TestCase):

    def test_job_initialization(self):
        job = Job()

if __name__ == '__main__':
   unittest.main()
