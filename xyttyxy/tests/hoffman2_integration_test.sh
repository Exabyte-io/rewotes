#!/bin/bash -f
#$ -cwd
#$ -o $JOB_ID.log
#$ -e $JOB_ID.err
#$ -pe dc* 16
#$ -l h_data=4G,h_rt=24:00:00

source /u/local/Modules/default/init/modules.sh # Source
source ~/.jobrc
module purge
module load intel/2020.4 intel/mpi

export VASP_PP_PATH=$HOME/bin/vasp_5.4.4/POT/
export OMP_NUM_THREAD=1
export I_MPI_COMPATIBILITY=4
export VASP_BIN='/u/home/x/xyttyxy/selfcompiled-programs/compile/vasp-5.4.4/source/bin/vasp_std_i2020.4_impi_vsol'
export VASP_COMMAND="mpirun -np 16 $VASP_BIN"

semiconductors=(Si_primitive Si_conventional GaAs CuO)
insulators=(SiO2 Al2O3 Fe2O3)
# SiO2 must decrease eps to 1e-3, otherwise will not converge
metals=(Cu Ni Au)

for material in ${semiconductors[@]}; do
    echo $material
done
		
for material in ${insulators[@]}; do
    echo $material
done
		
for material in ${metals[@]}; do
    echo $material
done



