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
dE=1 # (eV)



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
# Read the Total energy from the output


# Loop testing for dE

## set up calculation
#fileout='test.scf.in'
#qeinput.save(filein,fileout)

## run pw.x 

## extract the total energy

## compare with the threshold



