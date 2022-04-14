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
from pandas import DataFrame
from .data.bsl import pre_defined_libraries, nwchem_supported
from .qc import run_qc
import numpy as np
from sklearn.metrics import mean_absolute_percentage_error
from io import StringIO
import os
from .mol_classes import Atom,Molecule

def mol_from_xyz(xyz_file):
    mol = Molecule()
    with open(xyz_file,'r') as f:
        lines = f.readlines()
    for i in range(2,len(lines)):
        l = lines[i]
        if len(l.strip())>0:
            mol.add_atom(Atom(l.split()))

class BasisSetData():

    def __init__(self,name):
        self.name = name
        if self.name in nwchem_supported:
            self.supported_atom_types=nwchem_supported[name]
        else:
            self.supported_atom_types = None
        self.results = {}
        self.error = None
    
    def add_result(self,mol,res):
        self.results[mol.id] = res

    def calc_error(self,ref_data):
        pred = []
        ref = []
        for key in self.results:
            pred+=self.results[key]
            ref+=ref_data[key]
        self.error = mean_absolute_percentage_error(ref,pred)


class BasisSetOptimizer:

    def __init__(self,basis_library='double-zeta',prop_type='homo lumo gap',tolerance=0.01,functional='b3lyp'):
        self._ref_df = {}
        self._molecules = {}
        self._basis_library = self._setup_basis_library(basis_library)
        assert prop_type in ['homo lumo gap']
        self._prop_type = prop_type
        self._tolerance = tolerance
        self._functional = functional

    def _setup_basis_library(self,basis_library):
        _basis_library = {}
        if basis_library in pre_defined_libraries.keys():
            for basis in pre_defined_libraries[basis_library]:
                _basis_library[basis] = BasisSetData(basis)
            return _basis_library
        else:
            try:
                iter(basis_library)
            except TypeError:
                raise RuntimeError("Basis library should be one of the pre-defined libraries, or a list of basis sets.")
            for basis in basis_library:
                if basis in nwchem_supported:
                    _basis_library[basis] = BasisSetData(basis)
                else:
                    try:
                        with open(basis) as f:
                            f.next()
                        _basis_library[basis] = BasisSetData(basis)
                    except FileNotFoundError:
                        raise RuntimeError("{} is neither a supported basis set nor a path to a file.".format(basis))
            return _basis_library

    def _add_molecule(self,mol,ref_data):
        mol_id = 0
        while mol_id in self._molecules:
            mol_id+=1
        mol.id = mol_id
        self._molecules[mol_id] = mol
        self._ref_df[mol_id] = ref_data

    def add_molecule(self,xyz,ref_data,mol_units='Ang',prop_units='eV'):
        #TODO handle units
        if os.path.isfile(xyz):
            mol = mol_from_xyz(xyz)
        elif type(xyz)==str:
            mol = mol_from_xyz(StringIO(xyz))
        else:
            raise RuntimeError("Invalid specification of molecule. Must be a path to an .xyz file, or a string in .xyz format.")
        self._add_molecule(mol,ref_data)

    def optimize(self,functional=None,tolerance=None,engine='nwchem'):
        if tolerance is not None:
            self._tolerance = tolerance
        if functional is not None:
            self._functional = functional
        failed = {}
        success = []
        for name,basis in self._basis_library.items():
            for id,mol in self._molecules.items():
                res = run_qc(basis,mol,self._functional,self._prop_type,engine)
                if res['success']:
                    basis.add_result(mol,res['result'])
                else:
                    failed[name] = mol
                    break
            success.append(basis)
        for basis in success:
            basis.calc_error(self._ref_df)
        success = [x for x in success if x.error<self._tolerance]
        return {basis.name:basis.error for basis in sorted(success,key=lambda x: x.error)}
    

