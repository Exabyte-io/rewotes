import unittest
import andrewsalij.io as io
import os
from pymatgen.io import pwscf
import andrewsalij.Tests.test_params as test_params
QE_JOB = io.QEJob(test_params.FILE_PATH)
class test_qcinput_io(unittest.TestCase):
    def test_k_points(self):
        k_points = [1,2,3]
        QE_JOB.update_k_points(k_points)
        self.assertTupleEqual(tuple(k_points),QE_JOB.input.kpoints_grid)
    def test_write(self):
        try:
            os.makedirs(os.sep.join((test_params.TEST_DIR,"save")),exist_ok=True)
            QE_JOB.save(os.sep.join((test_params.TEST_DIR,"save/si2.in")))
        except:
            self.fail("Exception raised in saving QEJob")
    def test_run(self):
        try:
            QE_JOB.run()
        except:
            self.fail("Exception raised in running QEJob")
    def test_load_qejob(self):
        try:
            io.load_job(test_params.FILE_PATH,job_type="pwscf")
        except:
            self.fail("Exception raised in loading QEJob")
    def test_output(self):
        QE_JOB.run()
        assert isinstance(QE_JOB.output,pwscf.PWOutput)
if __name__ == '__main__':
    unittest.main()