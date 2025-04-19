#!/bin/bash

module load python
conda activate tdep
source /pscratch/sd/v/vladygin/doped-Si_project/MLFF_TDEP/testbench/bin/activate

export SLURM_CPU_BIND="cores"
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

module load PrgEnv-nvidia
module load cuda
module unload cudatoolkit/12.2
#module load nvhpc/23.1
module load cray-libsci/23.02.1.1
module load cray-hdf5-parallel
module swap gpu cpu 

python ../ConvergenceTracker/driver.py
