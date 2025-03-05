# Basis set selector

This package finds a basis set which can compute a chosen property with a given precision on a set of molecules.

# Requirements

In addition to a working installation of NWChem, this package requires the python libraries ase, rdkit, and numpy.


# Usage
First, use the script "select_basis_set.py" to select a basis set (or several) which meets your criteria. Usage of this script is as follows:


    python select_basis_set.py --input_files_path=test --property=energy --tolerance=0.1 --reference_values_file=example_energies.txt


If you do not provide an argument to "reference_values_file", it will compute the reference values usind CCSD.

When complete, this script will print out a list of basis sets that meets your tolerance requirement. 
These basis sets can be used as an argument to "calculate_property.py":

    python calculate_property.py --input_file=example.xyz --property=energy --basis_set='3-21G'

