from calculation import Calculation, VaspCalculation
from utils import ConvergenceProperty, messages
from error_metric import ErrorMetricScalar
from ase.io import read
from ase.calculators.vasp import Vasp
import os
import pytest


@pytest.fixture
def si_primitive_path():
    testdata = os.getenv("REWOTE_TEST_DATA")
    return os.path.join(testdata, "si_primitive")


@pytest.fixture
def si_primitive_atoms(si_primitive_path):
    return read(os.path.join(si_primitive_path, "POSCAR"))


@pytest.fixture
def si_primitive_calculator(si_primitive_path):
    return Vasp(directory=si_primitive_path, restart=True)


def test_calculation_abstract(si_primitive_atoms, si_primitive_path):
    with pytest.raises(TypeError):
        Calculation(
            "test_name",
            ConvergenceProperty.etotal,
            si_primitive_atoms,
            si_primitive_path,
        )


def test_vasp_calculation(si_primitive_atoms, si_primitive_path):
    try:
        VaspCalculation(
            "test_name",
            ConvergenceProperty.etotal,
            si_primitive_atoms,
            si_primitive_path,
        )
    except TypeError as ex:
        assert (
            False
        ), f"VaspCalculation constructor raised an exception {ex} when path {calc_path} is passed in"


def test_vasp_calculation_ref(si_primitive_atoms, si_primitive_calculator):
    try:
        VaspCalculation(
            "test_name",
            ConvergenceProperty.etotal,
            si_primitive_atoms,
            calculator=si_primitive_calculator,
        )

    except TypeError as ex:
        assert (
            False
        ), f"VaspCalculation constructor raised an exception {ex} when calculator is passed in"


@pytest.fixture
def si_primitive_calculation(si_primitive_atoms, si_primitive_path):
    calculation = VaspCalculation(
        "test_name", ConvergenceProperty.etotal, si_primitive_atoms, si_primitive_path
    )
    return calculation


def test_vasp_calculation_kpoints(si_primitive_calculation):
    # ensure kpoints info not read in
    # actually kpoints should probably be read by default and then cleared in KpointConvergenceTracker
    assert si_primitive_calculation.kpoints == (1, 1, 1)


def test_vasp_calculation_no_tetrahedron_smearing(si_primitive_calculation):
    calc = si_primitive_calculation._calculator
    assert calc.int_params["ismear"] == -5
    si_primitive_calculation.no_tetrahedron_smearing()
    assert calc.int_params["ismear"] == 0


def test_vasp_calculation_no_kspacing(si_primitive_calculation):
    calc = si_primitive_calculation._calculator
    assert calc.float_params["kspacing"] == 0.022
    si_primitive_calculation.no_kspacing()
    assert calc.float_params["kspacing"] == None


def test_vasp_calculation_set_to_singlepoint(si_primitive_calculation):
    calc = si_primitive_calculation._calculator
    assert calc.int_params["ibrion"] == 2
    assert calc.int_params["nsw"] == 500
    si_primitive_calculation.set_to_singlepoint()
    assert calc.int_params["ibrion"] == -1
    assert calc.int_params["nsw"] == 0


def test_vasp_calculation_calculation_required(si_primitive_calculation):
    with pytest.raises(Exception, match=messages("no_etotal")):
        e = si_primitive_calculation.etotal


@pytest.fixture
def si_primitive_calculation_222(si_primitive_atoms, si_primitive_path):
    # bypassed the restriction to not set etotal directly
    path = os.path.join(si_primitive_path, "2_2_2")
    calculator = Vasp(directory=path, restart=True)
    calculation = VaspCalculation(
        path,  # avoids resetting calculator
        ConvergenceProperty.etotal,
        si_primitive_atoms,
        calculator=calculator,
    )
    calculation.calculation_required = False
    return calculation


@pytest.fixture
def si_primitive_calculation_444(si_primitive_atoms, si_primitive_path):
    path = os.path.join(si_primitive_path, "4_4_4")
    calculator = Vasp(directory=path, restart=True)
    calculation = VaspCalculation(
        path, ConvergenceProperty.etotal, si_primitive_atoms, calculator=calculator
    )
    calculation.calculation_required = False
    return calculation


def test_error_metric_scalar(
    si_primitive_calculation_222, si_primitive_calculation_444
):
    error_metric = ErrorMetricScalar(fractional=False)
    error = error_metric.error(
        si_primitive_calculation_222, si_primitive_calculation_444
    )

    assert abs(error - 0.15795956) < 1e-10


def test_error_metric_scalar_fractional(
    si_primitive_calculation_222, si_primitive_calculation_444
):
    error_metric = ErrorMetricScalar(fractional=True)
    error = error_metric.error(
        si_primitive_calculation_222, si_primitive_calculation_444
    )
    assert abs(error - 0.014776226780188476) < 1e-10
