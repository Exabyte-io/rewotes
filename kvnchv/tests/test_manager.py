"""Preliminary pytest for components, to be broken up."""
import shutil
from pathlib import Path

import pytest
from converger import Manager
from converger.exceptions import SolverNotImplementedError, UnsupportedParameterError

resources = Path(__file__).parent.joinpath("resources")


input_v = {
    "solver_input": {
        "name": "vasp",
    },
    "input_path": str(resources),
    "tol": 1.0e-10,
    "target": "etot",
    "parameter_space": []
}

two_param_space_dict = {
    "solver_input": {
        "name": "espresso",
        "solver_path": "/usr/local/bin/pw.x"
    },
    "input_path": str('/home/chuk/testrun/'),
    "target": "etot",
    "tol": 1.0e-10,
    "parameter_space": [
        {
            "name": "k",
            "start": 1,
            "max": 5,
            "delta": 1
        },
        {
            "name": "p",
            "start": 0,
            "max": 10,
            "delta": 2
        },
    ]
}


@pytest.fixture
def two_param_job():
    """Job with unsupported parameters."""
    return Manager(two_param_space_dict)


# Test exceptions
def test_solver_not_implemented():
    """Test exception raised on 'vasp' solver name."""
    with pytest.raises(SolverNotImplementedError):
        Manager(input_v).run()


def test_two_params_unsupported(two_param_job):
    """Test unsupported parameter exception raised through manager."""
    with pytest.raises(UnsupportedParameterError):
        two_param_job.run()


def test_two_parameter_space(two_param_job):
    """Test generated parameter pairs for two parameters."""
    two_param_job.process_parameters()
    assert len(two_param_job.data) == 30


@pytest.fixture
def job_si(tmp_path):
    """Return valid espresso input object."""
    # copy test resources to tmp_path to avoid pollution
    shutil.copytree(resources.joinpath("Si"), tmp_path.joinpath("Si"))
    input_espresso_si = {
        "solver_input": {
            "name": "espresso",
        },
        "input_path": str(tmp_path.joinpath("Si")),
        "target": "etot",
        "tol": 1.0e-2,
        "parameter_space": [
            {
                "name": "k",
                "start": 1,
                "max": 5,
                "delta": 1
            },
        ]
    }
    return Manager(input_espresso_si)


def test_valid_input(job_si):
    """Test schema validation."""
    assert job_si.validate_input()


def test_solver_run(job_si):
    """Test Si job does not error and converges in 3 steps."""
    pytest.assume(job_si.run() == 0)
    pytest.assume(len(job_si.data.dropna()) == 3)
