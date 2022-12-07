# This code will run a k-point convergence up to a desire delta energy


import ioqeclass as qe

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

# maximum of iterations
Nitermax=20
# increasing step of kgrid
kstep=2

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
# Run initial test
##### pw.x test.scf.in test.scf.out
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

    # Run initial test
    ##### pw.x test.scf{counter}.in test.scf{counter}.out

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


