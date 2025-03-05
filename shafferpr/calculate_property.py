import argparse
import os, sys
import numpy as np
from BasisSetSelector import Molecule, Calculator

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Read in a molecular structure from an xyz file and compute a property of the molecule with a specified basis set ')
    parser.add_argument("--input_file", dest='input_file', default='.', help='Name of the xyz file containg the molecular structure.')
    parser.add_argument("--property", dest='property', help='The property you wish to calculate, select from "energy", "vibrational_energy".')
    parser.add_argument("--basis_set", dest='basis_set', help='The basis set you wish to use')
    args=parser.parse_args()

    molecule=Molecule(input_file_path=args.input_file,label=input_file)

    calculator=Calculator(basis_set=basis_set, property=args.property)
    value=calculator.calculate_property(molecule)
    print("The value of property %s for this molecule is %s"%(args.property,value))