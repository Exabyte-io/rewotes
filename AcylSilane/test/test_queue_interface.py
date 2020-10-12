import os
import unittest
from unittest import mock

from AcylSilane.convtrack.queue import Queue



class qstat_tests(unittest.TestCase):
    # Classwide qstat, for mocking subprocess.check_output('qstat')
    qstat_path = os.path.join(os.path.dirname(os.getcwd()), "data/sample_qstat.txt")
    with open(qstat_path, "r") as qstat:
        sample_qstat = qstat.read()

    def setUp(self):
        self.queue = Queue()

    def tearDown(self):
        pass

    @mock.patch("AcylSilane.convtrack.queue.subprocess.check_output",
                return_value="")
    def test_qstat_returns_empty_list_if_no_jobs(self, mock_subprocess):
        self.assertEqual(self.queue.qstat, [])

    @mock.patch("AcylSilane.convtrack.queue.subprocess.check_output",
                return_value=sample_qstat)
    def test_qstat_correct_length(self, mock_subprocess):
        expected_len = 2
        self.assertEqual(len(self.queue.qstat), 2)

    @mock.patch("AcylSilane.convtrack.queue.subprocess.check_output",
                return_value=sample_qstat)
    def test_qstat_parsed(self, mock_subprocess):
        expected_vals = {"Job-ID": "76174.master-production-20160630-cluster-001.exabyte.io",
                         "Username": "jrd101",
                         "Queue": "OR",
                         "Jobname": "VASP-TEST",
                         "SessID": "--",
                         "State": "Q",
                         "Nodes": "1:ppn=36",
                         "Req'd Time": "00:10:00",
                         "Rem'g Time": "--"}
        for key, val in expected_vals.items():
            for job in self.queue.qstat:
                self.assertEqual(job[key],val)

if __name__ == "__main__":
    unittest.main()
