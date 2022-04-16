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
from .data.bsl import pre_defined_libraries, nwchem_supported
from .qc import run_qc
from sklearn.metrics import mean_absolute_percentage_error
from io import StringIO
import os
from .mol_classes import Atom, Molecule
from typing import Union
import logging
log = logging.getLogger(__name__)

def mol_from_xyz(xyz_lines: str) -> Molecule:
    ''' Reads molecule from xyz formatted string.

    Parameters
    -----------
    xyz_lines: str
        XYZ specification of molecule.

    Returns
    --------
    mol: `~BasisSetSelector.mol_classes.Molecule`
        Molecule object

    '''
    mol = Molecule()
    for i in range(2, len(xyz_lines)):
        l = xyz_lines[i]
        if len(l.strip()) > 0:
            name, x, y, z = l.split()
            mol.add_atom(Atom(name, x, y, z))
    return mol


class BasisSetData():
    ''' Manages the computed data for a given basis set.

    Parameters
    -----------
    name: str
        Name of basis set
    '''

    def __init__(self, name: str):
        self.name = name
        self.results = {}
        self.error = None

    def add_result(self, mol_id: str, res: list):
        ''' Adds a computed result for a molecule.

        Parameters
        -----------
        mol_id: str
            ID of molecule
        res: list of float
            Result of calculation

        '''
        self.results[mol_id] = res

    def calc_error(self, ref_data: dict):
        ''' Computes mean absolute percent error of computed data relative to reference data.

        Parameters
        ----------
        ref_data: dict of str:list of float
            Reference data for each molecule
        '''

        pred = []
        ref = []
        for key in self.results:
            pred += self.results[key]
            ref += ref_data[key]
        self.error = mean_absolute_percentage_error(ref, pred)


class BasisSetOptimizer:
    ''' Determines the best basis set for a quantum chemistry calculation based on reference data for a chosen property.

    The user can specify which basis sets to consider, or choose one of the pre-curated lists. 

    Parameters
    -----------
    basis_library: str or list of str
        One of the pre-curated lists: 'double-zeta', 'triple-zeta', or a list of basis sets to consider.
    prop_type: str
        Property to optimize basis set for. Options: 'homo lumo gap', 'frequencies'.
    '''

    def __init__(self,
                 basis_library: Union[str, list] = 'double-zeta',
                 prop_type: str = 'homo lumo gap') -> None:
        self._ref_df = {}
        self._molecules = {}
        self._basis_library = self._setup_basis_library(basis_library)
        assert prop_type in ['homo lumo gap', 'frequencies']
        self._prop_type = prop_type

    def _setup_basis_library(self, basis_library:Union[str, list]):
        ''' Sets up basis set database.

        Parameters
        -----------
        basis_library: str or list of str
            One of the pre-curated lists: `double zeta`, `triple zeta`, or a list of basis sets to consider.

        '''
        _basis_library = {}
        if type(basis_library)==str:
            if basis_library in pre_defined_libraries.keys():
                for basis in pre_defined_libraries[basis_library]:
                    _basis_library[basis] = BasisSetData(basis)
                return _basis_library
        else:
            for basis in basis_library:
                _basis_library[basis] = BasisSetData(basis)
            return _basis_library

    def _add_molecule(self, mol: Molecule, ref_data: list) -> None:
        ''' Adds molecule and reference data point to the dataset.
        Parameters
        -----------
        mol: Molecule
            Molecule to add
        ref_data: list of float
            1D list of reference data points
        '''
        mol_id = 0
        while mol_id in self._molecules:
            mol_id += 1
        mol.id = mol_id
        self._molecules[mol_id] = mol
        self._ref_df[mol_id] = ref_data

    def add_basis_set(self,basis_set:str)->None:
        ''' Adds basis set to the basis library.
        Parameters
        -----------
        basis_set: str
            Name of basis set to add
        '''
        self._basis_library[basis_set] = BasisSetData(basis_set)

    def remove_basis_set(self,basis_set:str)->None:
        ''' Removes basis set from the basis library.
        Parameters
        -----------
        basis_set: str
            Name of basis set to remove
        '''
        self._basis_library.pop(basis_set, None)

    def add_molecule(self, xyz_file: str,
                              ref_data: Union[float, list]) -> None:
        ''' Adds molecule from .xyz file along with reference data point to optimizer.

        Parameters
        -----------
        xyz_file: str
            Path to .xyz file
        ref_data: float or list of float
            Reference data point(s). 1D data (e.g. frequencies) should be passed as a list of floats.
        '''
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        mol = mol_from_xyz(lines)
        if type(ref_data) == float:
            self._add_molecule(mol, [ref_data])
        else:
            self._add_molecule(mol, ref_data)

    def optimize(self,
                 functional: str,
                 precision: float,
                 engine: str = 'nwchem',
                 verbose:bool = False) -> dict:
        ''' Find basis sets which satisfy the target precision for the chosen functional.

        Parameters
        -----------
        precision: float
            Desired % error of computed data relative to reference data
        functional: str
            DFT functional name. Must be supported NWChem DFT functional
        engine: str
            Only NWChem supported for now.
        verbose: bool, optional
            Set to True to activate logging

        Returns
        --------
        result: dict of str:float
            Sorted dictionary of basis sets which satisfy the specified precision. 
        '''
        self._precision = precision
        self._functional = functional
        if verbose:
            logging.basicConfig(level=logging.INFO,format='%(message)s')
        failed = {}
        success = []
        for name, basis in self._basis_library.items():
            log.info('Running jobs for basis: {}'.format(name))
            for id, mol in self._molecules.items():
                log.info('Running molecule: {}'.format(mol.name))
                res = run_qc(basis, mol, self._functional, self._prop_type,
                             engine)
                if res['success']:
                    basis.add_result(mol.id, res['result'])
                else:
                    log.info("Job failed.")
                    failed[name] = mol
                    break
            if name not in failed:
                basis.calc_error(self._ref_df)
                success.append(basis)
        success = [x for x in success if x.error < self._precision]
        return {
            basis.name: basis.error
            for basis in sorted(success, key=lambda x: x.error)
        }
