#!/bin/sh

mpirun -n $SLURM_NPROCS singularity exec --bind HPL.dat:/HPL.dat ./hpl1.sif /bin/xhpl

