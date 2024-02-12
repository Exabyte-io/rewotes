import unittest
import andrewsalij.convergence_tracker as convergence_tracker

class CreateConvergenceTester(unittest.TestCase):
    def test_base_creation(self):
        convergence_tester = convergence_tracker.ConvergenceTester(path="")
        assert isinstance(convergence_tester,convergence_tracker.ConvergenceTester)
    def test_k_pt_creation(self):
        k_point_convergence_tester = convergence_tracker.KPointConvergenceTester(path="")
        assert isinstance(k_point_convergence_tester,convergence_tracker.KPointConvergenceTester)
if __name__ == '__main__':
    unittest.main()
