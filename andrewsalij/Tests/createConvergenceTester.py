import unittest
import andrewsalij.convergence_tracker as convergence_tracker
import andrewsalij.io
import andrewsalij.Tests.test_params as test_params

class CreateConvergenceTester(unittest.TestCase):
    def test_base_creation(self):
        convergence_tester = convergence_tracker.ConvergenceTester(path="")
        assert isinstance(convergence_tester,convergence_tracker.ConvergenceTester)
    def test_k_pt_creation(self):
        k_point_convergence_tester = convergence_tracker.KPointConvergenceTester(path="")
        assert isinstance(k_point_convergence_tester,convergence_tracker.KPointConvergenceTester)
    def test_convergence_loaded(self):
        k_point_convergence_tester = convergence_tracker.KPointConvergenceTester(test_params.FILE_PATH)
        assert isinstance(k_point_convergence_tester.base_job,andrewsalij.io.Job)

if __name__ == '__main__':
    unittest.main()
