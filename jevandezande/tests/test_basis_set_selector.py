# isort: skip_file
from pytest import approx, fixture, mark

from jevandezande.basis_set_selector import energy, get_homo_lumo_gap, optimize_basis_set

H2_geom = "H\nH 1 1"
H2O_geom = "H\nO 1 1\nH 2 1 1 105"


@fixture
def H2():
    return energy(H2_geom)


@fixture
def H2O():
    return energy(H2O_geom)


@mark.parametrize(
    "mol, e_tot",
    [
        ["He", -2.89812055849662],
        [H2_geom, -1.150637241322931],
        [H2O_geom, -76.35782896158842],
    ],
)
def test_energy(mol, e_tot):
    assert energy(mol).e_tot == approx(e_tot)


def test_get_homo_lumo_gap(H2, H2O):
    assert get_homo_lumo_gap(H2) == approx(-0.33167556339619)
    assert get_homo_lumo_gap(H2O) == approx(-0.2528650709452739)


@mark.parametrize(
    "mol, target, tolerance, success, basis_set",
    [
        ["He", 1, 0.1, False, "def2-SVP"],
        ["He", -1.6416, 0.001, True, "def2-SVP"],
        [H2_geom, 1, -0.1, False, "def2-SVP"],
        [H2_geom, -0.3280, 0.001, True, "def2-TZVP"],
        [H2O_geom, -1, 0.000000001, False, "def2-SVP"],
        [H2O_geom, -0.2529, 0.001, True, "def2-SVP"],
    ],
)
def test_optimize_basis_set(mol, target, tolerance, success, basis_set):
    energy_kwargs = {
        "geom": mol,
        "functional": "BP86",
    }
    best = optimize_basis_set(
        energy_kwargs,
        get_homo_lumo_gap,
        target=target,
        tolerance=tolerance,
        stop_early=True,
    )

    print(best)
    assert best["success"] == success
    assert best.get("basis_set", None) == basis_set
