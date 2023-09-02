# Overview

> When passed a directory path containing the POSCAR, KPOINTS, INCAR and POTCAR VASP files, as well as a convergence threshold value, finds the number of KPOINT divisions that converge the DFT total energy within the convergence threshold
> Usage: ``` python run.py -d {DIRECTORY PATH} -t {CONVERGENCE_THRESHOLD} ```
> Outputs: Printed table showing the path to KPOINTS convegence, as well as the number of divisions that converge within the supplied threshold 

# Dependencies
> numpy
> ASE (Atomic Simulation Environment)
> tabulate

# Convergence
> Begins all calculations from [1, 1, 1]. Incrementally increases the number of divisions along each reciprocal lattice vector by decreasing the KSPACING between divisions. Keeps a constant KSPACING value for each calculation, as recommended by VaspWiki. 

# Example (Si2)

