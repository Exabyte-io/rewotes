import argparse
import os, sys
import numpy as np
from BasisSetSelector import Molecule, Calculator


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Read a set of molecules from xyz files, and select a basis set that allows you to calculate a chosen property with a desired tolerance')
    parser.add_argument("--input_files_path", dest='input_files_path', default='.', help='path of directory containing the xyz files')
    parser.add_argument("--property", dest='property', help='the property you wish to optimize, select from "energy", "vibrational_energy" or ')
    parser.add_argument("--tolerance", dest='tolerance', type=float, help='the tolerance you wish to maintain for calculating the property of interest, expressed as a percentage')
    parser.add_argument("--reference_values_file", dest='reference_values_file', help= 'a file with a list of reference values which will be use to assess the quality of each basis set,\
                        if this file is not provided, the reference values will be calculated automatically with CCSD')
    args=parser.parse_args()

    # read files from input directory
    input_files=os.listdir(args.input_files_path)
    molecules=[]
    for file in input_files:
        molecules.append(Molecule(input_file_path=args.input_files_path+'/'+file,label=file))


    basis_sets=['STO-2G','STO-3G','STO-6G','3-21G','4-31G','6-31G'] #list of basis sets

    #populate a dictionary with reference values either from a file, or by doing calculations using ccsd
    reference_values_dict={}
    if args.reference_values_file is not None:
        with open(args.reference_values_file,"r") as f:
            x=f.readlines()
        reference_values_dict={q.split()[0] : float(q.split()[1]) for q in x}
    else:
        calculator=Calculator(property=args.property, use_reference_calculator=True)
        reference_values_dict={x.label : calculator.calculate_property(x) for x in molecules}


    #loop over basis sets and calculate the property of interest for each one.
    accurate_basis_sets=[]
    for basis_set in basis_sets:
        calculator=Calculator(basis_set=basis_set,property=args.property)
        values_dict={x.label : calculator.calculate_property(x) for x in molecules}
        percent_deviation_dict = {key : np.abs((values_dict[key]-reference_values_dict[key])/reference_values_dict[key]) for key in values_dict}
        maximum_deviation=max(percent_deviation_dict.values())
        if maximum_deviation < args.tolerance:
            accurate_basis_sets.append(basis_set)
        #print(percent_deviation_dict)

    print("The following basis sets are sufficiently accurate:\n")
    print(accurate_basis_sets)