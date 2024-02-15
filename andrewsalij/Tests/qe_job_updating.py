import copy
import unittest
import test_params
import andrewsalij.io
import numpy as np
import numpy.testing as np_test
QE_JOB = andrewsalij.io.QEJob(test_params.FILE_PATH)
class QEJobUpdateTester(unittest.TestCase):
    def test_k_array_update(self):
        new_k_points = np.array([1,7,4])
        QE_JOB.update_k_points(new_k_points)
        np_test.assert_array_equal(new_k_points,np.array(QE_JOB.input.kpoints_grid))
    def test_section_update(self):
        new_directory = "/home/new_directory"
        section_dict_update = {'control':{'outdir':new_directory}}
        QE_JOB_copy = copy.deepcopy(QE_JOB)
        QE_JOB_copy.update_input_sections(section_dict_update)
        new_outdir = QE_JOB_copy.input.sections['control']['outdir']
        assert new_outdir == new_directory
if __name__ == '__main__':
    unittest.main()
