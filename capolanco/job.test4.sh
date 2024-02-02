#!/bin/bash

#PBS -N test
#PBS -l nodes=1
#PBS -l ppn=8
#PBS -q OR
#PBS -j oe
#PBS -l walltime=05:03:00

cd $PBS_O_WORKDIR
module load espresso
mpirun -np $PBS_NP pw.x test.scf4.in test.scf4.out
