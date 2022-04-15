#!/bin/sh

mpirun -n $SLURM_NNODES singularity exec "$@"
#--bind HPL.dat:/HPL.dat ./hpl1.sif /bin/xhpl

