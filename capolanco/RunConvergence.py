# This code will run a k-point convergence up to a desire delta energy

import subprocess
import ioqeclass as qe
import ioclusterclass as cluster

###########################
##### USER INPUTS
###########################

# Define the path of the input file
# The k-grid of this file will be set as the starting point
# for the convergence process
filein='Si.scf.in'

# Define delta energy threshold for 
# convergence (eV)
dEthreshold=1.0 # (eV)


###########################
##### DEVELOPERS INPUTS
###########################

## Convergence
# maximum of iterations
Nitermax=20
# increasing step of kgrid
kstep=2

## Cluster
# number of nodes to be used
Nnodes=1
# number of processors per node
ppn=8
# queue
queue='OR'
# walltime
walltimehours=5
walltimeminutes=3


###########################
##### PROGRAM
###########################


### Set up initial parameters
# Initialize pw.x input class
qeinput=qe.qepwinput()
# load the pw.x input 
qeinput.load(filein)
# Create the input file for initial run
testin='test.scf.in'
testout='test.scf.out'
qeinput.save(filein,testin)
# Create the job to send to the cluster
# Initialize job class
job=cluster.jobclass()
# set up job class
job.name='test'
job.nodes=Nnodes
job.ppn=ppn
job.queue=queue
job.walltimehours=walltimehours
job.walltimeminutes=walltimeminutes
# create the job file
jobname=f'job.test.sh'
job.createjobQEpw(jobname,testin,testout)
# Run initial test
subprocess.run(['echo',f'runing {jobname}'])
##subprocess.run(['qsub',f'{jobname}'])

# Initialize pw.x output class
qeoutput=qe.qepwoutput()
# Read the Total energy from the output
qeoutput.getenergy(testout)


# Loop testing for dE
dE=2.0*dEthreshold
EnergyOld=qeoutput.energy
counter=0
while ((dE>dEthreshold) and (counter<Nitermax)):

    # Increase counter
    counter=counter+1
    print(f'\n## Iteration {counter}')

    # Increase k grid
    qeinput.kgrid=qeinput.kgrid+kstep

    # Create the input file for initial run
    testin=f'test.scf{counter}.in'
    testout=f'test.scf{counter}.out'
    qeinput.save(filein,testin)

    # create the job file
    jobname=f'job.test{counter}.sh'
    job.createjobQEpw(jobname,testin,testout)
    # Run QE calculation
    subprocess.run(['echo',f'runing {jobname}'])
    ##subprocess.run(['qsub',f'{jobname}'])

    # Read the Total energy from the output
    qeoutput.getenergy(testout)

    # Update dE and EnergyOld
    dE=abs(EnergyOld-qeoutput.energy)
    EnergyOld=qeoutput.energy
    print(f'dE {dE} eV')

# Display results
if (dE<dEthreshold):
    print(f'Convergence has been achieved in {counter} iterations')
    print(f'for kgrid {qeinput.kgrid}')
    print(f'The total enegy change is less than {dEthreshold} eV')
else:
    print(f'Convergence has NOT been achieved in {counter} iterations')


