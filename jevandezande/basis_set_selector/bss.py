from collections.abc import Collection

import pyscf

AHLRICHS_BASIS_SETS = (
    "def2-SVP",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-TZVPD",
    "def2-TZVPPD",
    "def2-QZVP",
    "def2-QZVPD",
    "def2-QZVPPD",
)


def energy(
    geom: str,
    functional: str = "BP86",
    basis_set: str = "def2-SVP",
    name: str = "molecule",
):
    mol = pyscf.M(atom=geom, basis=basis_set, symmetry=False)

    mf = mol.RKS()
    mf.xc = functional
    mf.run()

    return mf


def get_homo_lumo_gap(mf, thresh: float = 1e-10) -> float:
    for i, occ in enumerate(mf.mo_occ):
        if occ < thresh:
            HOMO = i - 1
            LUMO = i
            break

    return mf.mo_energy[HOMO] - mf.mo_energy[LUMO]


def optimize_basis_set(
    energy_kwargs: dict,
    property_function,
    target: float,
    tolerance: float = 0.01,
    basis_sets: Collection = AHLRICHS_BASIS_SETS,
    stop_early: bool = True,
) -> dict:
    """
    Find a basis set within the targetted relative error.

    :param energy_kwargs: parameters for the energy function. `basis_set` will be modified for each run.
    :param property_function: 1D property function to run
    :param target: value to optimize for
    :param tolerance: maximum relative_error allowed
    :param basis_sets: basis sets to optimize with
    :param stop_early: stop testing once a satisficing basis set has been found
    """
    best = {"success": False}

    return best


if __name__ == "__main__":
    H2 = "H 0 0 0\nH 0 0 1"

    mf = energy(H2)
    hl_gap = get_homo_lumo_gap(mf)
    print(hl_gap)
