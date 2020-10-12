class Calculation(object):
    def __init__(self, incar, poscar, kpoints=None):
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints