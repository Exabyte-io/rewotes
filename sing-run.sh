#!/bin/sh

mpirun -n $SLURM_NPROCS singularity exec "$@"

