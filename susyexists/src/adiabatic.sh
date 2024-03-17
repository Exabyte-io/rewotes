#!/bin/zsh

#### Local testing variables ####
#Path to quantum-espresso
QE=/home/susy/q-e/bin/
#Number of cores
SLURM_NTASKS=4

#Initial degauss value (0 for POSCAR initialization)
initial=0
#POSCAR file
POSCAR=C.poscar
project_name=Graphene
pw_json=pw.json
#vc-relax k points
vc_relax_k=('1 1 1 0 0 0')
relax_k=('1 1 1 0 0 0')
scf_k=('1 1 1 0 0 0')


echo "Espresso Machine is starting"
#Smearing parameters
for sigma in 0.05 
do
echo "Adiabatic automation for $sigma is started"

#Crystal relaxation
if [ $initial -eq 0 ]
then 
    echo "Generating crystal relaxation input from POSCAR for $sigma"
    python run.py -n $project_name -c vc-relax -j $pw_json  -d $sigma -p $POSCAR -k ${vc_relax_k[@]}
else 
    echo "Generating crystal relaxation input from output for $sigma"
    python run.py -n $project_name -c vc-relax -j $pw_json  -d $sigma -i ./$project_name/$initial/vc-relax.out -k ${vc_relax_k[@]}
fi

echo "Calculating crystal relaxation for $sigma"
mpirun -np $SLURM_NTASKS $QE/pw.x -inp ./$project_name/$sigma/vc-relax.in > ./$project_name/$sigma/vc-relax.out
echo "Crystal relaxation for $sigma is finished"

#Atom Relaxation
echo "Generating relaxation input for $sigma"
python run.py -n $project_name -c relax -j $pw_json  -d $sigma -k ${relax_k[@]} -l mono
echo "Starting relaxation for $sigma"
mpirun -n $SLURM_NTASKS $QE/pw.x -inp ./$project_name/$sigma/relax.in > ./$project_name/$sigma/relax.out
echo "Atom relaxation for $sigma is finished"

#Scf calculation
echo "Scf input for $sigma is generated"
python run.py -n $project_name -c scf -j $pw_json  -d $sigma -k ${scf_k[@]}
echo "Scf calculation for $sigma is started"
mpirun -n $SLURM_NTASKS $QE/pw.x -inp ./$project_name/$sigma/scf.in > ./$project_name/$sigma/scf.out
echo "Scf calculation for $sigma is finished"

#Band calculation
echo "Band calculation input for $sigma is generated"
python run.py -n $project_name -c bands -j $pw_json  -d $sigma
echo "Band calculation for $sigma is started"
mpirun -n $SLURM_NTASKS $QE/pw.x -inp ./$project_name/$sigma/bands.in > ./$project_name/$sigma/bands.out
echo "Band-pp for $sigma is generated"
python run.py -n $project_name -c bands-pp -j $pw_json  -d $sigma
echo "Band-pp calculation $sigma is started"
$QE/bands.x -inp ./$project_name/$sigma/bands-pp.in > ./$project_name/$sigma/bands-pp.out
mkdir ./$project_name/plots/ 2> /dev/null
python ./src/plot_band.py -n $project_name  -d $sigma

echo "Adiabatic automation for $sigma is completed"
##Memorize the last calculation
initial=$sigma
done 

echo "Adiabatic automation is completed"
echo "vc-relax > relax > scf > band"

# # mkdir dyn
# # mkdir dyn/$sigma/
# # python run.py -n $project_name -c ph -j $pw_json  -d $sigma
# # mpirun -n 256 $QE/ph.x -npool 256 -inp ./inputs/$sigma-ph.in > ./outputs/$sigma-ph.out
# # python run.py -n $project_name -c q2r -j $pw_json  -d $sigma
# # $QE/q2r.x  < ./inputs/$sigma-q2r.in > ./outputs/$sigma-q2r.out
# # python run.py -n $project_name -c matdyn -j $pw_json  -d $sigma
# # $QE/matdyn.x  < ./inputs/$sigma-matdyn.in > ./outputs/$sigma-matdyn.out
# # python run.py -n $project_name -c plotband -j $pw_json  -d $sigma
# # cd ./dyn/$sigma/
# # $QE/plotband.x  < $sigma-plotband.in > $sigma-plotband.out
# # cd ../..
# # bash ph.sh $QE $sigma
# python ph_plot.py -j $pw_json  -d $sigma
