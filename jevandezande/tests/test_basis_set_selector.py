from pytest import approx, fixture, mark

from jevandezande.basis_set_selector import get_homo_lumo_gap, run


@fixture
def H2():
    return run("H 0 0 0\nH 0 0 1")


@fixture
def H2O():
    return run("H\nO 1 1\nH 2 1 1 105")


@mark.parametrize(
    "mol, energy",
    [
        ["He", -2.89812055849662],
        ["H 0 0 0\nH 0 0 1", -1.150637241322931],
        ["H 0 0 0\nO 0 0 1\nH 0 1 1", -76.35496969473999],
    ],
)
def test_run(mol, energy):
    assert run(mol).e_tot == approx(energy)


def test_get_homo_lumo_gap(H2, H2O):
    assert get_homo_lumo_gap(H2) == approx(-0.33167556339619)
    assert get_homo_lumo_gap(H2O) == approx(-0.2528650709452739)
