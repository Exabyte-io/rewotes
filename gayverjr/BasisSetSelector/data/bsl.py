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
# Gaussian basis sets pulled from: https://nwchemgit.github.io/AvailableBasisSets.html
nwchem_supported = {
    '3-21++G': [
        'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '3-21++G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    '3-21G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs'
    ],
    '3-21G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    '3-21GSP': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    '4-22GSP': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    '4-31G':
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'P', 'S', 'Cl'],
    '6-31++G': [
        'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-31++G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-31++G**': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-31+G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-311++G(2d,2p)': [
        'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-311++G(3df,3pd)': [
        'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar'
    ],
    '6-311++G**': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-311+G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    '6-311G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br',
        'Kr', 'I'
    ],
    '6-311G(2df,2pd)':
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'K', 'Ca'],
    '6-311G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br',
        'Kr', 'I'
    ],
    '6-311G**': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br',
        'Kr', 'I'
    ],
    '6-31G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    '6-31G-Blaudeau': ['K', 'Ca'],
    '6-31G(3df,3pd)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    '6-31G*': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    '6-31G*-Blaudeau': ['K', 'Ca'],
    '6-31G**': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'Ahlrichs pVDZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'Ahlrichs TZV': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'Ahlrichs VDZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'Ahlrichs VTZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'ANO-RCC': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'
    ],
    'apr-cc-pV(Q+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-cc-pCV5Z':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pCVDZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar'
    ],
    'aug-cc-pCVQZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar'
    ],
    'aug-cc-pCV(T+d)Z': ['Al', 'Si', 'P', 'S', 'Cl'],
    'aug-cc-pCVTZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar'
    ],
    'aug-cc-pV(5+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pV5Z': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
        'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-cc-pV(6+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pV6Z': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-cc-pV(D+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pVDZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-cc-pV(Q+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pVQZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-cc-pV(T+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pVTZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-cc-pVTZ-J': [
        'H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl', 'Sc', 'Ti',
        'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Se'
    ],
    'aug-cc-pwCV5Z':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pwCV5Z-NR':
    ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'aug-cc-pwCVDZ':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pwCVQZ':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Br'],
    'aug-cc-pwCVQZ-NR':
    ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'aug-cc-pwCVTZ':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'aug-cc-pwCVTZ-NR':
    ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'aug-mcc-pV5Z': ['H'],
    'aug-mcc-pV6Z': ['H'],
    'aug-mcc-pV7Z': ['H'],
    'aug-mcc-pV8Z': ['H'],
    'aug-mcc-pVQZ': ['H'],
    'aug-mcc-pVTZ': ['H'],
    'aug-pc-0': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br',
        'Kr'
    ],
    'aug-pc-1': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-pc-2': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-pc-3': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-pc-4': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'aug-pcJ-0': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-0_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-1': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-1_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-2': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-2_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-3': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-3_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-4': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcJ-4_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'aug-pcS-0': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-pcS-1': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-pcS-2': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-pcS-3': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-pcS-4': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'aug-pV7Z': ['C', 'N', 'O', 'F', 'S'],
    'B2 basis set for Zn': ['Zn'],
    'Bauschlicher ANO': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    'Binning/Curtiss SV': ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    'Binning/Curtiss SVP': ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    'Binning/Curtiss VTZ': ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    'Binning/Curtiss VTZP': ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    'cc-pCV5Z': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar', 'Ca'
    ],
    'cc-pCV6Z':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pCV6Z(old)': ['O'],
    'cc-pCVDZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'Ca'
    ],
    'cc-pCVDZ(old)': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pCVQZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'Ca'
    ],
    'cc-pCVQZ(old)': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pCVTZ': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'Ca'
    ],
    'cc-pCVTZ(old)': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pV(5+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pV5Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pV(6+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pV6Z': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'cc-pV8Z': ['H', 'Ne'],
    'cc-pV9Z': ['Ne'],
    'cc-pV(D+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pVDZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pVDZ(seg-opt)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pV(Q+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pVQZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pVQZ(seg-opt)': ['H', 'He', 'C', 'O'],
    'cc-pV(T+d)Z': ['Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pVTZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pVTZ(seg-opt)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'cc-pwCV5Z':
    ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'cc-pwCV5Z Core Set': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'cc-pwCV5Z-NR': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'cc-pwCVDZ': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'cc-pwCVQZ': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Br'
    ],
    'cc-pwCVQZ-NR': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'cc-pwCVTZ': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'cc-pwCVTZ-NR': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'ccemd-2': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'ccemd-3': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'ccJ-pV5Z': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'ccJ-pVDZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'ccJ-pVQZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'ccJ-pVTZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'coemd-2': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'coemd-3': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'coemd-4': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'coemd-ref': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'CVTZ': ['Li', 'Be', 'Na', 'Mg', 'K', 'Ca'],
    'd-aug-cc-pV5Z': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'd-aug-cc-pV6Z': ['H', 'B', 'C', 'N', 'O'],
    'd-aug-cc-pVDZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'd-aug-cc-pVQZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'd-aug-cc-pVTZ': ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'Def2-QZVP': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Def2-QZVPD': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Def2-SVP': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Def2-SVPD': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Def2-TZVP': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Def2-TZVPD': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'dhf-QZVP': [
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'dhf-SV(P)': [
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'dhf-TZVP': [
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
    ],
    'Dunning-Hay Double Rydberg': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'Dunning-Hay Rydberg': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'DZ + Double Rydberg (Dunning-Hay)': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl'
    ],
    'DZ + Rydberg (Dunning)': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl'
    ],
    'DZ (Dunning)': [
        'H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl'
    ],
    'DZP + Rydberg (Dunning)': [
        'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl'
    ],
    'DZP (Dunning)': [
        'H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl'
    ],
    'DZQ': ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag'],
    'DZVP2 (DFT Orbital)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S',
        'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'DZVP (DFT Orbital)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe'
    ],
    'Feller Misc. CVDZ': ['K'],
    'Feller Misc. CVQZ': ['K'],
    'Feller Misc. CVTZ': ['K'],
    'GAMESS PVTZ': [
        'H', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar'
    ],
    'GAMESS VTZ': [
        'H', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar'
    ],
    'IGLO-II': ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl'],
    'IGLO-III': ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl'],
    'jul-cc-pV(D+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'jul-cc-pV(Q+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'jul-cc-pV(T+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'jun-cc-pV(D+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'jun-cc-pV(Q+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'jun-cc-pV(T+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'LANL08': [
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
        'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
        'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf',
        'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi'
    ],
    'LANL08+': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'LANL08d': [
        'Si', 'P', 'S', 'Cl', 'Ge', 'As', 'Se', 'Br', 'Sn', 'Sb', 'Te', 'I',
        'Pb', 'Bi'
    ],
    'LANL08(f)': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb',
        'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os',
        'Ir', 'Pt', 'Au'
    ],
    'Lanl2-[10s8p7d3f2g]': ['Rh'],
    'Lanl2-[5s4p4d2f]': ['Rh'],
    'Lanl2-[6s4p4d2f]': ['Rh'],
    'Lanl2DZ+1d1f': ['Rh'],
    'Lanl2DZ+2s2p2d2f': ['Rh'],
    'LANL2TZ': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr',
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'La', 'Hf', 'Ta', 'W',
        'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
    ],
    'LANL2TZ+': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'LANL2TZ(f)': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb',
        'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os',
        'Ir', 'Pt', 'Au'
    ],
    'm6-31G': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    'm6-31G*': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    'maug-cc-pV(D+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'maug-cc-pVDZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'maug-cc-pV(Q+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'maug-cc-pVQZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'maug-cc-pV(T+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'maug-cc-pVTZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'may-cc-pV(Q+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'may-cc-pV(T+d)Z': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'McLean/Chandler VTZ': ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'MG3S': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'MIDI!': ['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'],
    'MIDI (Huzinaga)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K', 'Cs'
    ],
    'MINI (Huzinaga)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    'MINI (Scaled)': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    ],
    'modified LANL2DZ': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb',
        'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os',
        'Ir', 'Pt', 'Au'
    ],
    'NASA Ames ANO': [
        'H', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'P', 'Ti', 'Fe', 'Ni'
    ],
    'NASA Ames ANO2': ['Sc', 'Ti', 'V'],
    'NASA Ames cc-pCV5Z': ['Ti'],
    'NASA Ames cc-pCVQZ': ['Ti'],
    'NASA Ames cc-pCVTZ': ['Ti'],
    'NASA Ames cc-pV5Z': ['Ti', 'Fe'],
    'NASA Ames cc-pVQZ': ['Ti', 'Fe'],
    'NASA Ames cc-pVTZ': ['Ti', 'Fe'],
    'Partridge Uncontracted 1': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr'
    ],
    'Partridge Uncontracted 2': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'Partridge Uncontracted 3': [
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
        'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn'
    ],
    'Partridge Uncontracted 4': [
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'pc-0': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br',
        'Kr'
    ],
    'pc-1': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'pc-2': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'pc-3': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'pc-4': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'pcemd-2': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'pcemd-3': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'pcemd-4': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'pcJ-0': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-0_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-1': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-1_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-2': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-2_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-3': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-3_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-4': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcJ-4_2006': [
        'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl',
        'Ar'
    ],
    'pcS-0': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'pcS-1': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'pcS-2': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'pcS-3': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'pcS-4': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'pSBKJC': [
        'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Ge', 'As', 'Se', 'Br', 'Sn',
        'Sb', 'Te', 'I'
    ],
    'Pt - mDZP': ['Pt'],
    'pV6Z': ['H', 'C', 'N', 'O'],
    'pV7Z': ['H', 'C', 'N', 'O', 'F', 'Ne', 'S'],
    'Roos_ANO_DZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'Roos_ANO_TZ': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn'
    ],
    'Roos Augmented Double Zeta ANO': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn'
    ],
    'Roos Augmented Triple Zeta ANO': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
        'Ni', 'Cu', 'Zn'
    ],
    's3-21G': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    's3-21G*': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    's6-31G': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    's6-31G*': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'Sadlej pVTZ': [
        'H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Si', 'P', 'S', 'Cl',
        'K', 'Ca', 'Br', 'Rb', 'Sr', 'I'
    ],
    'SDB-aug-cc-pVQZ': [
        'Ga', 'Ge', 'As', 'Se', 'Br', 'In', 'Sn', 'Sb', 'Te', 'I'
    ],
    'SDB-aug-cc-pVTZ': [
        'Ga', 'Ge', 'As', 'Se', 'Br', 'In', 'Sn', 'Sb', 'Te', 'I'
    ],
    'SDB-cc-pVQZ': [
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'
    ],
    'SDB-cc-pVTZ': [
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'
    ],
    'STO-2G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sr'
    ],
    'STO-3G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I'
    ],
    'STO-3G*': ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    'STO-6G': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
    ],
    'SV + Double Rydberg (Dunning-Hay)': ['B', 'C', 'N', 'O', 'F', 'Ne'],
    'SV + Rydberg (Dunning-Hay)': ['B', 'C', 'N', 'O', 'F', 'Ne'],
    'SV (Dunning-Hay)': ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'SVP + Rydberg (Dunning-Hay)': ['B', 'C', 'N', 'O', 'F', 'Ne'],
    'SVP (Dunning-Hay)': ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'TZ (Dunning)': ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
    'TZVP (DFT Orbital)': [
        'H', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    ],
    'UGBS': [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pu', 'Am', 'Cf', 'Es', 'Fm', 'Md', 'No',
        'Lr'
    ],
    'un-ccemd-ref': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'un-pcemd-ref': [
        'He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar',
        'P', 'Cl', 'Al', 'Si', 'S'
    ],
    'Wachters+f': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    'WTBS': [
        'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn'
    ],
    'Z3Pol': ['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl']
}

pre_defined_libraries = { 'double-zeta' : ['6-31G','6-31+G*','6-31++G**','cc-pVDZ','aug-cc-pVDZ',
                                           'Def2-SVP','Def2-SVPD','pc-2','aug-pc-2'],
                          'triple-zeta' : ['6-311G','6-311+G*','6-311++G**','cc-pVTZ','aug-cc-pVTZ',
                                           'Def2-TZVPD','Def2-TZVPPD','pc-3','aug-pc-3']}