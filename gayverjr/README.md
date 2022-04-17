# Basis set selector (Chemistry)

> Ideal candidate: scientists skilled in Density Functional Theory and proficient in python.

# Overview

The aim of this task is to create a simple python package that implements automatic basis set selection mechanism for a quantum chemistry engine.

# Requirements

1. automatically find the basis set delivering a particular precision, passed as argument (eg. within 0.01% from reference)
1. use either experimental data or higher-fidelity modeling results (eg. coupled cluster) as reference data
1. example properties to converge: HOMO-LUMO gaps, vibrational frequencies

# Expectations

- mine reference data for use during the project
- correctly find a basis set that satisfies a desired tolerance for a set of 10-100 molecules, starting from H2, as simplest, up to a 10-20-atom ones
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# User story

As a user of this software I can start it passing:

- molecular structure
- reference datapoint
- tolerance (precision)

as parameters and get the basis set that satisfies the tolerance criterion.

# Notes

- create an account at exabyte.io and use it for the calculation purposes
- suggested modeling engine: NWCHEM or SIESTA


# Discussion

- I decided to implement the package as both a command line script (optimize_basis) and a python package (see examples directory).
- JSON files seemed like a straightforward way to address the challenge of specifying geometries/data for a set of 10-100 molecules 
- In addition to allowing users to specify a list of basis sets, the package also allows users to choose pre-curated lists of basis sets (double-zeta and triple-zeta).
- I chose to assume that the user has NWChem installed on their machine and used temporary directories to manage the scratch files. The `qc.py` file can in principle be modified to work with other engines.
- The reference data and geometries were obtained from the CCCBDB, and was CCSD/aug-cc-pVQZ unless otherwise specified.
