# K-point convergence tracker (Materials)

Dependencies for this project may be installed by running in the following install script.

```bash
pip install numpy matplotlib 
pip install pymatgen
pip install pydantic 
```

Tests (see Scripts folder) have been run on a variety of materials, whose provenance is below:

Materials Project
link: https://legacy.materialsproject.org/
(DOI: 10.1063/1.4812323)
(CC-BY 4.0, https://creativecommons.org/licenses/by/4.0/legalcode), see https://legacy.materialsproject.org/terms

Si2 (mp-149) (DOI 10.17188/1190959) (https://legacy.materialsproject.org/materials/mp-149/)

Materials Cloud three-dimensional crystals database (MC3D) 
link: https://archive.materialscloud.org/record/2022.38
(DOI: 10.24435/materialscloud:rw-t0)
(CC-BY 4.0, https://creativecommons.org/licenses/by/4.0/legalcode)


GaN (mc3d-3763/pbe)
BN (mc3d-13290/pbe)
O4Ru2 (mc3d-1930/pbe)
Br4Ca2 (mc3d-30836/pbe)
Cs2La2Te6Zn2 (mc3d-11071/pbe)
F2Fe2O8S2 (mc3d-14425/pbe)

Materials Cloud two-dimensional crystals database (MC2D)
link: https://archive.materialscloud.org/record/2020.158
DOI:10.24435/materialscloud:az-b2 
DOI:10.24435/materialscloud:36-nd (expansion)
(CC-BY 4.0, https://creativecommons.org/licenses/by/4.0/legalcode)

C (graphene, from graphite exfoliation) (https://www.materialscloud.org/discover/mc2d/details/C, graphite (2H) initial)
MoS2 (https://www.materialscloud.org/discover/mc2d/details/MoS2-MoS2)
AgCO2 (https://www.materialscloud.org/discover/mc2d/details/AgCO2)


Pseudopotentials Tested:
PBEsol (standard accuracy) NC SR ONCVPSP 0.4.1 from Pseudo Dojo 
	Paper: 10.1016/j.cpc.2018.01.012 arxiv preprint 
	Method Paper: 10.1103/PhysRevB.88.085117
	License: (CC-BY 4.0, https://creativecommons.org/licenses/by/4.0/legalcode) 
	(see https://github.com/PseudoDojo/pseudodojo)
		

> Ideal candidate: scientists skilled in Density Functional Theory and proficient in python.

# Overview

The aim of this task is to create a python package that implements automatic convergence tracking mechanism for a materials simulations engine. The convergence is tracked with respect to the k-point sampling inside a reciprocal cell of a crystalline compound.

# Requirements

1. automatically find the dimensions of a k-point mesh that satisfy a certain criteria for total energy (eg. total energy is converged within dE = 0.01meV)
1. the code shall be written in a way that can facilitate easy addition of convergence wrt other characteristics extracted from simulations (forces, pressures, phonon frequencies etc)
1. the code shall support VASP or Quantum ESPRESSO

# Expectations

- correctly find k-point mesh that satisfies total energy convergence parameters for a set of 10 materials, starting from Si2, as simplest, to a 10-20-atom supercell of your choice
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# User story

As a user of this software I can start it passing:

- path to input data (eg. pw.in / POSCAR, INCAR, KPOINTS) and
- kinetic energy cutoff

as parameters and get the k-point dimensions (eg. 5 5 5).

# Notes

- create an account at exabyte.io and use it for the calculation purposes
- suggested modeling engine: Quantum ESPRESSO
