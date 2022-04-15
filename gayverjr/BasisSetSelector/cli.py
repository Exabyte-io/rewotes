import argparse
import json
from .mol_classes import Atom,Molecule
from .Optimizer import BasisSetOptimizer
from typing import Tuple
import logging
parser = argparse.ArgumentParser(description='Optimize basis set for a set a of molecules+reference data for a given property.')
parser.add_argument('json_file', metavar='JSON File', type=str,help='JSON file specifying job information.')
parser.add_argument(
    '-v', '--verbose',
    help="Verbose output.",
    action="store_const", dest="loglevel", const=logging.INFO,
    default=logging.WARNING,
)

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
    functional: str
        Chosen DFT functional
    tolerance: float
        Chosen % error threshold
    '''
    with open(json_file, "r") as f:
        data = json.load(f)
    sections = data['Basis Set Selector']
    molecules = []
    basis_library = None
    ref_data = None
    prop_type = None
    tolerance = None
    functional = None
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
            tolerance = section['tolerance']
            functional = section['functional']
        else:
            raise RuntimeError("Not a valid section.")
    opt = BasisSetOptimizer(basis_library,prop_type)
    for i,mol in enumerate(molecules):
        if type(ref_data[i]) == float:
            opt._add_molecule(mol, [ref_data[i]])
        else:
            opt._add_molecule(mol, ref_data[i])
    return opt, functional, tolerance

def main():
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel,format='%(message)s')
    opt, tolerance, functional = read_json(args.json_file)
    result = opt.optimize(tolerance,functional)
    if result == {}:
        print("None of the selected basis sets satisfy the chosen tolerance.")
    else:
        print("Results:")
        for basis,error in result.items():
            print('{}: {}'.format(basis,error))