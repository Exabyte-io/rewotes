convergence test
================

'convergence test' rewote for Willis O'Leary (wolearyc)

Automatically determines minimal k-point sampling to achieve a well-converged property (e.g. energy) 

usage: python convergencetracker.py [-h] [-i K K K] [-c CORES] [-i ECUTS] dir threshold
-i defines the intial k-point grid from which to begin the convergence test
-c defines the number of cores to use in calculations
-i defines the desired kinetic energy cutoffs in Ry for the wavenumbers and charge density/potential
dir specifies a directory which must contain an appropriate pw.in file
	You can specify k-points in this file, though they will be changed over the course of the convergence test.
	outdir, wtcdir, and pseudo_dir must be specified with templates (see the Si example)
	The ATOMIC_POSITIONS, CELL_PARAMETERS, and K_POINTS blocks must follow immediately after each other in that order (no spaces). 
threshold defines the energy convergence threshold in eV

When run, the program will first initialize a new material in Mat3ra based on the structural information found in the pw.in file. It will then run several jobs with progressively finer k-point grids. In the current implementation, only k-point grids with identical sampling ratios as the initial k-point grid are evaluated. For example, with an initial grid of 1x2x2, the program will evaluate 2x4x4, 3x6x6, etc.

If one does not specify an initial k-point grid, the default starting k-point grid is a 1x1x1 grid (gamma-centered). Not only is this a rather poor starting guess for small supercells, but it will not evenly sample reciprocal space for many systems. Therefore, one is advised to carefully choose a reasonable initial k-point grid based on previous knowledge of the system (i.e. the literature) as well as the lengths of the (reciprocal) lattice vectors.

Currently, the program does not support shifted k-point grids.
