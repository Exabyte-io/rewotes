"""Single solver instance tests."""
import shutil
from pathlib import Path

import numpy as np
import pytest
from converger import EspressoSolver
from converger.exceptions import EspressoOutputNotFoundError, UnsupportedParameterError

resources = Path(__file__).parent.joinpath("resources")

supported_param = {"k": 1}

expected_Si_k1 = {'etot': -7.84373769449927, 'fnorm': 2.0744045342550352e-07}

unsupported_params = {"k": 1, "p": 2}

base_solver_input = {
        "name": "espresso",
        "solver_path": "/home/runner/q-e-qe-7.1/bin/pw.x"
    }


@pytest.fixture
def test_dir_tmp_si(tmp_path):
    """Return temporary testing directory for Si."""
    tmp_outdir = tmp_path.joinpath("Si")
    shutil.copytree(resources.joinpath("Si"), tmp_outdir)
    return tmp_outdir


def test_raises_unsupported_parameter_error(test_dir_tmp_si):
    """Test unsupported parameter raised directly."""
    with pytest.raises(UnsupportedParameterError):
        EspressoSolver(base_solver_input, test_dir_tmp_si, unsupported_params)


def test_solver_result_close_k1(test_dir_tmp_si):
    """Test EspressoSolver result for k=1."""
    solver_k1 = EspressoSolver(base_solver_input, test_dir_tmp_si, supported_param)
    solver_k1.run()
    np.testing.assert_allclose(list(solver_k1.results.values()),
                               list(expected_Si_k1.values()))


def test_raises_output_not_found_error(test_dir_tmp_si):
    """Test espresso output not found rasied when no output data."""
    solver_k1 = EspressoSolver(base_solver_input, test_dir_tmp_si, supported_param)
    with pytest.raises(EspressoOutputNotFoundError):
        solver_k1.parse_output(test_dir_tmp_si)
