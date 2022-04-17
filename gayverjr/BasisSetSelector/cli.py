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
import argparse
import json
from .mol_classes import Atom,Molecule
from .Optimizer import BasisSetOptimizer
from typing import Tuple
import logging
parser = argparse.ArgumentParser(description='Optimize basis set for a set a of molecules+reference data for a given property.')
parser.add_argument('json_file', metavar='JSON File', type=str,help='JSON file specifying job information.')
parser.add_argument('functional',metavar='DFT Functional', type=str, help='DFT functional to optimize basis set for.')
parser.add_argument('precision', metavar='Precision', type=float, help='Desired precision of basis set (percent error).')
parser.add_argument(
    '-v', '--verbose',
    help="Verbose output.",
    action="store_const", dest="loglevel", const=logging.INFO,
    default=logging.WARNING,
)

#TODO input validation for reading json
def read_json(json_file:str)->Tuple[BasisSetOptimizer,str,float]:
    ''' Parses JSON file to set up basis set optimization.

    Parameters
    -----------
    json_file: str
        JSON file 

    Returns
    --------
    opt: `~BasisSetSelector.BasisSetOptimizer`
        Optimizer for basis set
    '''
    with open(json_file, "r") as f:
        data = json.load(f)
    sections = data['Basis Set Selector']
    molecules = []
    basis_library = None
    ref_data = None
    prop_type = None
    for section in sections:
        if section['title'] == 'molecule':
            for mol in section['geometry']:
                new_mol = Molecule()
                for atom in mol:
                    atm_name = atom['atom']
                    xyz = atom['xyz']
                    new_mol.add_atom(Atom(atm_name,xyz[0],xyz[1],xyz[2]))
                molecules.append(new_mol)
        elif section['title'] == 'basis':
            basis_library = section['basis_library']
        elif section['title'] == 'reference':
            ref_data = section['data']
            prop_type = section['property']
        else:
            raise RuntimeError("Not a valid section.")
    opt = BasisSetOptimizer(basis_library,prop_type)
    for i,mol in enumerate(molecules):
        if type(ref_data[i]) == float:
            opt._add_molecule(mol, [ref_data[i]])
        else:
            opt._add_molecule(mol, ref_data[i])
    return opt

def main():
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel,format='%(message)s')
    opt = read_json(args.json_file)
    result = opt.optimize(args.functional,args.precision)
    if result == {}:
        print("None of the selected basis sets satisfy the chosen precision.")
    else:
        print("Results:")
        for basis,error in result.items():
            print('{}: {}'.format(basis,error))