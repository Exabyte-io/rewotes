"""Preliminary pytest for components, to be broken up."""
from pathlib import Path

from converger import Manager

import pytest

resources = Path(__file__).parent.joinpath("resources")

input_espresso = {
    "solver_input": {
        "name": "espresso",
        "solver_path": "/usr/local/bin/pw.x"
    },
    "input_path": str(resources),
    "tol": 1.0e-10,
    "target": "total_energy"
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
def manager():
    """Espresso input object."""
    return Manager(input_espresso)


def test_valid_input(manager):
    """Test schema validation."""
    assert manager.validate_input()


def test_solver_not_implemented():
    """Test exception raised on 'vasp' solver name."""
    with pytest.raises(Exception):
        Manager(input_v)


def test_solver_run(manager):
    """Test single solver run through manager class."""
    manager.run_solver()
