from pytest import approx, fixture, mark

from jevandezande.basis_set_selector import get_homo_lumo_gap, optimize_basis_set, run


H2_geom = "H\nH 1 1"
H2O_geom = "H\nO 1 1\nH 2 1 1 105"


@fixture
def H2():
    return run(H2_geom)


@fixture
def H2O():
    return run(H2O_geom)


@mark.parametrize(
    "mol, energy",
    [
        ["He", -2.89812055849662],
        [H2_geom, -1.150637241322931],
        [H2O_geom, -76.35496969473999],
    ],
)
def test_run(mol, energy):
    assert run(mol).e_tot == approx(energy)


def test_get_homo_lumo_gap(H2, H2O):
    assert get_homo_lumo_gap(H2) == approx(-0.33167556339619)
    assert get_homo_lumo_gap(H2O) == approx(-0.2528650709452739)


@mark.xfail
@mark.parametrize(
    "mol, target, tolerance, success",
    [
        ["He", 1, 0.1, False],
        ["He", 1, 1, True],
        [H2_geom, 1, 0.1, False],
        [H2_geom, 1.2, 0.1, True],
        [H2O_geom, 1, 0.01, False],
    ],
)
def test_optimize_basis_set(mol, target, tolerance, success):
    run_kwargs = {
        "geom": mol,
        "functional": "BP86",
    }
    best = optimize_basis_set(
        run_kwargs,
        get_homo_lumo_gap,
        target=target,
        tolerance=tolerance,
        stop_early=True,
    )
    assert best.success == success
