import argparse
import os, sys

from BasisSetSelector import Molecule


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Read a set of molecules from xyz files, and select a basis set that allows you to calculate a chosen property with a desired tolerance')
    parser.add_argument("--input_files_path", dest='input_files_path', default='.', help='path of directory containing the xyz files')
    parser.add_argument("--property", dest='property', help='the property you wish to optimize, select from "energy", "vibrational_energy" or ')
    parser.add_argument("--tolerance", dest='tolerance', help='the tolerance you wish to maintain for calculating the property of interest, expressed as a percentage')
    parser.add_argument("--reference_values_file", dest='reference_values_file', help= 'a file with a list of reference values which will be use to assess the quality of each basis set,\
                        if this file is not provided, the reference values will be calculated automatically with CCSD')
    args=parser.parse_args()

    input_files=os.listdir(args.input_files_path)
    molecules=[]
    for file in input_files:
        molecules.append(Molecule(input_file_path=args.input_files_path+'/'+file))
    basis_sets=['3-21G','4-31G','6-31G']
    for basis_set in basis_sets:
        values=[x.calculate_energy(basis_set) for x in molecules]
        print(values)