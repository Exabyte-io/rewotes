from .data.elements import ELEMENTS
import numpy as np

class Atom(object):
    def __init__(self, name, x, y, z):
        assert name in ELEMENTS
        self.name = name
        self.charge = ELEMENTS[self.name]
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.xyz = np.array([self.x, self.y, self.z])

    def __str__(self):
        return '  {:}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format(self.name[0], self.xyz[0], self.xyz[1], self.xyz[2])

class Molecule(object):
    def __init__(self):
        self.atoms = []
        self.id = None

    def add_atom(self,atom):
        self.atoms.append(atom) 

    def __str__(self):
        string_output = ''
        for atom in self.atoms:
            string_output += atom.__str__()
        return string_output.rstrip()