#!/bin/bash -f
#$ -cwd
#$ -o $JOB_ID.log
#$ -e $JOB_ID.err
#$ -pe dc* 16
#$ -l h_data=4G,h_rt=24:00:00

source /u/local/Modules/default/init/modules.sh # Source
module add intel/17.0.7

source ~/.bashrc

export VASP_PP_PATH=$HOME/Potentials_5.4.1/POT_GGA_PAW_PBE/
export OMP_NUM_THREAD=1
export I_MPI_COMPATIBILITY=4
export VASP_COMMAND='mpirun -np ${NSLOTS} /u/project/sautet/shared/VASP.5.4.1/vasp.5.4.1/bin/vasp_std'
python geo-opt.py
echo "run complete on `hostname`: `date` `pwd`" >> ~/job.log
