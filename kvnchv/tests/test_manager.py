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
    "target": "total_energy",
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
    "target": "total_energy"
}


@pytest.fixture
def job():
    """Espresso input object."""
    return Manager(input_espresso)


def test_valid_input(job):
    """Test schema validation."""
    assert job.validate_input()


def test_solver_not_implemented():
    """Test exception raised on 'vasp' solver name."""
    with pytest.raises(Exception):
        Manager(input_v)


def test_solver_run(job):
    """Test single solver run through manager class."""
    job.run()


def test_two_parameter_space():
    """Test generated parameter pairs for two parameters."""
    two_param_space_dict = {
        "solver_input": {
            "name": "espresso",
            "solver_path": "/usr/local/bin/pw.x"
        },
        "input_path": str('/home/chuk/testrun/'),
        "target": "total_energy",
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
    job = Manager(two_param_space_dict)
    job._process_parameters()
    assert len(job.data) == 20
