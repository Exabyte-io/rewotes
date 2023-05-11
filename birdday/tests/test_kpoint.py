import pytest
from kpoint.kpoint import ConvTracker

def test_check_convergence():
    # Below convergence limit
    test_tracker = ConvTracker("placeholder_config", "placeholder_job_endpoint", cutoff=1e-5, energy=[5, 5+1e-6])
    assert test_tracker.check_convergence() == True

    # At convergence limit
    test_tracker = ConvTracker("placeholder_config", "placeholder_job_endpoint", cutoff=1e-5, energy=[5, 5+1e-5])
    assert test_tracker.check_convergence() == True

    # Above convergence limit
    test_tracker = ConvTracker("placeholder_config", "placeholder_job_endpoint", cutoff=1e-5, energy=[5, 5+1e-4])
    assert test_tracker.check_convergence() == False
