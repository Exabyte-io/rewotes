# Overview

* When passed a directory path containing the POSCAR, KPOINTS, INCAR and POTCAR VASP files, as well as a convergence threshold value, finds the number of KPOINT divisions that converge the DFT total energy within the convergence threshold
* Usage: ``` python run.py -d {DIRECTORY PATH} -t {CONVERGENCE_THRESHOLD} ```
* Outputs: Printed table showing the path to KPOINTS convegence, as well as the number of divisions that converge within the supplied threshold 

# Dependencies
* numpy for kpoints division calculations
* ASE (Atomic Simulation Environment) for VASP job submission
* tabulate for printing the KPOINT convergence trajectory

# Convergence
* Begins all convergence calculations at [1, 1, 1].
* Finds the largest KSPACING value that increases the number of kpoint divisions along a reciprocal lattice vector by 1: $$max([KSPACING_{1}, KSPACING_{2}, KSPACING_{3}])$$
  and calculates the number of divisions according to this value.
* Iteratively decreases the KSPACING value used for each calculation regardless of the reciprocal lattice geometry, as recommended by VaspWiki, until the convergence threshold is reached. 

# Example (Si2)
![Screen Shot 2023-09-05 at 5 09 27 PM](https://github.com/rymo1354/rewotes/assets/52838869/f4ca1bbf-64f7-48c5-a14e-5a880a539630)


