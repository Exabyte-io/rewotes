import argparse
import json
from .mol_classes import Atom,Molecule
from .Optimizer import BasisSetOptimizer
from pandas import DataFrame
parser = argparse.ArgumentParser(description='Parse a JSON file.')
parser.add_argument('json_file', metavar='JSON File', type=str,help='JSON file specifying job information.')

def read_json(json_file):
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
    opt = BasisSetOptimizer(basis_library,prop_type,tolerance,functional)
    for i,mol in enumerate(molecules):
        opt._add_molecule(mol,ref_data[i])
    return opt

def main():
    args = parser.parse_args()
    opt = read_json(args.json_file)
    result = opt.optimize()
    if result == {}:
        print("None of the selected basis sets satisfy the chosen tolerance.")
    else:
        print("Results:")
        for basis,error in result.items():
            print('{}: {}'.format(basis,error))