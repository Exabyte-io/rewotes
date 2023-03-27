"""Preliminary pytest for components, to be broken up."""
from pathlib import Path

import pytest
from converger import Manager


resources = Path(__file__).parent.joinpath("resources")

input_espresso = {
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
    ]
}

input_v = {
    "solver_input": {
        "name": "vasp",
    },
    "input_path": str(resources),
    "tol": 1.0e-10,
    "target": "etot"
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
def job():
    """Espresso input object."""
    return Manager(input_espresso)


@pytest.fixture
def two_param_job():
    """Job with unsupported parameters."""
    return Manager(two_param_space_dict)


def test_valid_input(job):
    """Test schema validation."""
    assert job.validate_input()


def test_solver_not_implemented():
    """Test exception raised on 'vasp' solver name."""
    with pytest.raises(Exception):
        Manager(input_v)


def test_two_params_unsupported(two_param_job):
    """Test exception raised on unsupported parameter."""
    with pytest.raises(Exception):
        two_param_job.run()


def test_solver_run(job):
    """Test single solver run through manager class."""
    job.run()


def test_two_parameter_space(two_param_job):
    """Test generated parameter pairs for two parameters."""
    two_param_job._process_parameters()
    assert len(two_param_job.data) == 20
