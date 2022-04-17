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

import subprocess
import os
from BasisSetSelector import BasisSetOptimizer
from BasisSetSelector import read_json
import numpy as np

dirname = os.path.dirname(__file__)
json_file = os.path.join(dirname,'data.json')
xyz_file = os.path.join(dirname,'h2.xyz')

def test_cli():
    output = subprocess.run(['optimize_basis', json_file,'b3lyp','0.1'], capture_output=True,encoding='UTF-8')
    assert '3-21G: 0.0212' in output.stdout

def test_constructors():
    opt = BasisSetOptimizer('triple-zeta','homo lumo gap')
    assert all(x in opt._basis_library for x in ['6-311G','6-311+G*','6-311++G**','cc-pVTZ','aug-cc-pVTZ',
                                           'Def2-TZVPD','Def2-TZVPPD','pc-3','aug-pc-3'])
    opt = read_json(json_file)
    assert all(x in opt._basis_library for x in ["3-21G","STO-3G"])

def test_add_basis():
    opt = BasisSetOptimizer(["3-21G","STO-3G"],'homo lumo gap')
    opt.add_basis_set('cc-pvdz')
    assert all(x in opt._basis_library for x in ["3-21G","STO-3G",'cc-pvdz'])

def test_optimize_frequencies():
    opt = read_json(json_file)
    result = opt.optimize('b3lyp',0.1)
    assert np.round(result['3-21G'],decimals=3) == 0.021
    result2 = opt.optimize('pbe0',0.1) 
    assert np.round(result2['3-21G'],decimals=3) == 0.022

def test_optimize_homo_lumo():
    opt = BasisSetOptimizer(["3-21G","STO-3G"],'homo lumo gap')
    opt.add_molecule(xyz_file,17.4479)
    result = opt.optimize('b3lyp',0.2)
    assert np.round(result['3-21G'],decimals=3) == 0.145
    assert 'STO-3G' not in result