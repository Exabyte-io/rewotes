import logging
import math
from collections.abc import Collection

import pyscf

from .basis_sets import AHLRICHS_BASIS_SETS


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

    if "basis_set" in energy_kwargs:
        raise ValueError(
            "`basis_set` cannot be adjusted as it is automatically set by optimize_basis_set()."
        )

    for basis_set in basis_sets:
        logging.info(f"{basis_set=}")

        mf = energy(**energy_kwargs, basis_set=basis_set)
        prop = property_function(mf)

        relative_error = (prop - target) / target
        logging.info(f"{relative_error=}")

        if abs(relative_error) < best.get("relative_error", math.inf):
            best |= {
                "basis_set": basis_set,
                "property": prop,
                "relative_error": relative_error,
            }
            if abs(relative_error) < tolerance:
                best["success"] = True
                logging.info(f"Met {target=} with {relative_error=} and {basis_set=}")
                if stop_early:
                    return best

    if best.get("relative_error", math.inf) > target:
        logging.info("Unable to find an acceptable basis set in {basis_sets}")

    return best


if __name__ == "__main__":
    H2 = "H 0 0 0\nH 0 0 1"

    mf = energy(H2)
    hl_gap = get_homo_lumo_gap(mf)
    print(hl_gap)

    energy_kwargs = {
        "geom": H2,
        "functional": "BP86",
    }
    best = optimize_basis_set(
        energy_kwargs,
        get_homo_lumo_gap,
        target=-0.3238,
        tolerance=0.001,
        stop_early=True,
    )
    for key, value in best.items():
        print(f"{key}: {value}")
