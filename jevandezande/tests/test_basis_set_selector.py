from pytest import approx, fixture, mark

from jevandezande.basis_set_selector import get_homo_lumo_gap, run


@fixture
def H2():
    return run("H 0 0 0\nH 0 0 1")


def test_run():
    mf = run("H 0 0 0\nH 0 0 1")
    assert mf.e_tot == approx(-1.12761841620146)


@mark.xfail
def test_get_homo_lumo_gap(H2):
    get_homo_lumo_gap(H2)
    assert 0
