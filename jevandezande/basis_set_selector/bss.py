import pyscf


def run(geom: str, functional: str = "BP86", basis_set: str = "def2-SVP", name: str = "molecule"):
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


if __name__ == "__main__":
    mf = run("H 0 0 0\nH 0 0 1")
    print(mf)
    hl_gap = get_homo_lumo_gap(mf)
    print(hl_gap)
