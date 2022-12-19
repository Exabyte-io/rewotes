# Supercell convergence tracker (Materials)

> Ideal candidate: scientists skilled in Density Functional Theory and proficient in python.

# Overview

Some DFT codes lack reciprocal-space capabilities. For example, [until 2015][1] CP2K did not have the capability to utilize k-points, and certain features (such as symmetry-based reductions of k-points) still remain to be implemented.
A common solution to this problem is to generate a supercell to represent the material.
The aim of this task is to create a python package that implements an automatic convergence tracking mechanism for a materials simulations engine. The convergence is tracked with respect to the supercell size of a crystalline material.

This task was modified from the k-point convergence-tracker ReWoTe.

# Requirements

1. automatically find the dimensions of a supercell that satisfy a certain criteria for total energy (eg. total energy is converged within dE = 0.01meV)
1. the code shall be written in a way that can facilitate easy addition of convergence wrt other characteristics extracted from simulations (forces, pressures, phonon frequencies etc)
1. the code shall support VASP

# Expectations

- correctly find the supercell size that satisfies total energy convergence parameters for a Cu FCC unit cell (note: task originally stated several types of cell here, but has been modified based on our conversation, to minimize the test's computational needs)
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# User story

As a user of this software I can start it passing:

- path to input data (eg. pw.in / POSCAR, INCAR, KPOINTS) and
- kinetic energy cutoff

as parameters and get the supercell dimensions (eg. 5 5 5).

# Notes

- create an account at exabyte.io and use it for the calculation purposes
- modeling engine: VASP

[1]: <https://www.cp2k.org/faq:kpoints> "CP2K FAQ: Which features are working with k-point sampling in CP2K?"