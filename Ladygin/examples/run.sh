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

curr_path=$( pwd  )

for mat in Al Si MbB2 NaCl Ni GaAs Cu2O ZrO2 SiO2 CsPdBr3 Pb; do
        rm -r calcs/${mat}
        cp -r inputs/${mat} calcs
	ConvTrack -workdir "${curr_path}/calcs/${mat}" -mode "qe" -target "total_energy" -eps 0.01 -input "pw.in" -encut 40 -calc "par" -ncores 32 -nk 8 -k_range 2 80 4

done
