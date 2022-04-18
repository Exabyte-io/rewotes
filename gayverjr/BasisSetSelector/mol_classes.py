## Copyright (c) 2022, James Gayvert
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from .data.elements import ELEMENTS

class Atom(object):
    def __init__(self, name, x, y, z):
        assert name in ELEMENTS
        self.name = name
        self.charge = ELEMENTS[self.name]
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.xyz = [self.x, self.y, self.z]

    def __str__(self):
        return '  {:}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format(self.name[0], self.xyz[0], self.xyz[1], self.xyz[2])

class Molecule(object):
    def __init__(self):
        self.atoms = []
        self.id = None

    def add_atom(self,atom):
        self.atoms.append(atom) 

    def get_natoms(self):
        return len(self.atoms)

    def get_name(self):
        return ''.join(x.name for x in self.atoms)

    natoms = property(get_natoms)

    name = property(get_name)

    def __str__(self):
        string_output = ''
        for atom in self.atoms:
            string_output += atom.__str__()
        return string_output.rstrip()