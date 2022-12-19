import os
import unittest
from unittest import mock

from AcylSilane.convtrack.queueinterface import QueueInterface


class qstat_tests(unittest.TestCase):
    # Classwide qstat, for mocking subprocess.check_output('qstat')
    qstat_path = os.path.join(os.getcwd(), "data/simulated_qstat.txt")
    with open(qstat_path, "r") as qstat:
        sample_qstat = qstat.read()

    def setUp(self):
        self.queue = QueueInterface()

    def tearDown(self):
        pass

    @mock.patch("AcylSilane.convtrack.queueinterface.subprocess.check_output",
                return_value="")
    def test_qstat_returns_empty_list_if_no_jobs(self, mock_subprocess):
        self.assertEqual(self.queue.qstat(), [])

    @mock.patch("AcylSilane.convtrack.queueinterface.subprocess.check_output",
                return_value=sample_qstat)
    def test_qstat_correct_length(self, mock_subprocess):
        expected_len = 4
        self.assertEqual(len(self.queue.qstat()), expected_len)

    @mock.patch("AcylSilane.convtrack.queueinterface.subprocess.check_output",
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
            for job in self.queue.qstat()[0:2]:
                self.assertEqual(job[key], val)


class qsub_tests(unittest.TestCase):
    # Classwide qsub, for mocking subprocess.check_output('qstat')
    qsub_path = os.path.join(os.getcwd(), "data/simulated_qsub.txt")
    with open(qsub_path, "r") as qsub:
        sample_qsub = qsub.read()

    def setUp(self):
        self.queue = QueueInterface()

    def tearDown(self):
        pass

    @mock.patch("AcylSilane.convtrack.queueinterface.subprocess.check_output",
                return_value=sample_qsub)
    def test_qsub_returns_jobid(self, mock_qsub):
        job_id = self.queue.qsub("job_vasp.pbs")
        mock_qsub.assert_called_with(["qsub", "job_vasp.pbs"])
        self.assertEqual(job_id, self.sample_qsub)
