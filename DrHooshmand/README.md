# Generalized Stacking Fault Energy (GSFE) Surface Generator


# Overview

The aim of this task is to create a python package that implements automatic generation of GSFE in hexagonal close packed (HCP) materials. 

# User story

As a user of this software, I can predict the landscape of the energy penalty
between two adjacent planes during shear deformation in a specific slip
direction on a given slip plane in a target material. The result of this analysis represents the nature of slip and invloves the stable and unstable stacking fault energies. 


# Requirements

- Basic python packages such as numpy, math and os are necessary. 
- The code shall generate all displacement needed to generate the input files
  for capturing the GSFE landscape
- These input files will be used in molecular dynamics (MD) or density
  functional theory (DFT) codes to calculate the energy terms and construct the
  GSFE surface

# Input/Output

- The inputs required for the code: 
    1) Conventional supercell in which the displacements are going to be applied
    2) Range or displacement steps that are required to be applied to the cell
    3) Number of periodic units in each direction of fault plane (XY) 

- Output: Generated POSCAR (or LAMMPS/data) files corresponding to each
  displacement

# Expectations

- the code shall be able to suggest realistic values for GSFE
- test case would be provided


# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- use a designated github repository for version control
