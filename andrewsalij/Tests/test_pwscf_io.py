import unittest
import andrewsalij.io as io
import os
TEST_DIR = "/home/andrew/Documents/RewoteTests"
FILE_STR = "si2.in"
FILE_PATH = os.sep.join((TEST_DIR,FILE_STR))
QE_JOB = io.QEJob(FILE_PATH)
class test_qcinput_io(unittest.TestCase):
    def test_k_points(self):
        k_points = [1,2,3]
        QE_JOB.update_k_points(k_points)
        self.assertTupleEqual(tuple(k_points),QE_JOB.input.kpoints_grid)
    def test_write(self):
        try:
            os.makedirs(os.sep.join((TEST_DIR,"save")),exist_ok=True)
            QE_JOB.save(os.sep.join((TEST_DIR,"save/si2.in")))
        except:
            self.fail("Exception raised")
    def test_run(self):
        try:
            QE_JOB.run()
        except:
            self.fail("Exception raised")

if __name__ == '__main__':
    unittest.main()