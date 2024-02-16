import unittest
import numpy as np
from andrewsalij.convergence_tracker import KPointConvergenceTester

class KPointArrayCreation(unittest.TestCase):
    def test_basic(self):
        k_array = KPointConvergenceTester.create_k_point_array(k_iterator=1)
        assert np.array_equal(k_array,np.array([1,1,1],dtype=int))
    def test_k_indices(self):
        for i in np.arange(3):
            k_array = KPointConvergenceTester.create_k_point_array(k_iterator=2,k_index=i)
            assert np.array_equal(k_array,np.array([2,2,2],dtype=int))
    def test_weight(self):
        weight = np.array([1,2.4,3.1])
        k_array = KPointConvergenceTester.create_k_point_array(k_iterator=10,k_index = 0,effective_weight_array=weight)
        assert np.array_equal(k_array,np.array([10,24,31]))
    def test_k_force(self):
        weight = np.array([3,2.4,2])
        k_force_array = np.array([None,None,1]) #for 2D (xy) materials
        k_array = KPointConvergenceTester.create_k_point_array(k_iterator=5,k_index =1,effective_weight_array=weight,
                                                               fixed_k_points=k_force_array)
        assert np.array_equal(k_array, np.array([6, 5, 1]))
    def test_k_force_invalid_warn(self):
        k_force=  1 #should warn--scalar does not work
        self.assertWarns(Warning,KPointConvergenceTester.create_k_point_array,5,
                         k_index =1,effective_weight_array = np.array([3.6,3.1,6.3]),fixed_k_points = k_force)
    def test_k_force_invalid_pass(self):
        weight = np.array([3,2.4,2])
        k_force_array = 1 #should pass with warning
        k_array = KPointConvergenceTester.create_k_point_array(k_iterator=5,k_index =1,effective_weight_array=weight,
                                                               fixed_k_points=k_force_array)
        assert np.array_equal(k_array, np.array([6, 5, 4]))

if __name__ == '__main__':
    unittest.main()
