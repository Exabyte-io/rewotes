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
