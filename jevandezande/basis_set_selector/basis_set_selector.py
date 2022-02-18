from __future__ import annotations

from psi4 import core
from psi4.driver import energy, geometry

core.IOManager.shared_object()


def run(geom: str, method: str = "BP86", basis: str = "def2-SVP", name: str = "molecule"):
    geometry(geom + "\n\nsymmetry c1", name)

    core.IO.set_default_namespace(name)

    e, wfn = energy(f"{method}/{basis}", return_wfn=True)

    return wfn


def get_homo_lumo_gap(wfn) -> float | tuple[float, float]:
    a, b = wfn.nalpha(), wfn.nbeta()
    ea, eb = wfn.epsilon_a().to_array(), wfn.epsilon_b().to_array()

    if wfn.molecule().multiplicity() == 1:
        return ea[a + 1] - ea[a]

    return ea[a + 1] - ea[a], eb[b + 1] - eb[b]


if __name__ == "__main__":
    wfn = run("H\nH 1 1")
    hl_gap = get_homo_lumo_gap(wfn)
    print(hl_gap)
