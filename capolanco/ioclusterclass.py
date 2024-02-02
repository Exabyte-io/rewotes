import numpy as np


## THINGS TO FIX

class jobclass:
    """This class creates jobs for the cluster"""

    ################### 
    ## Initialization function for the class
    ###################
    def __init__(self):
        
        ## name of the job
        self.name='Unset'
        ## number of nodes
        self.nodes=np.int64(1)
        ## number of processors per node
        self.ppn=np.int64(1)
        ## queue
        self.queue='D'
        ## walltime (in hours)
        self.walltimehours=np.int64(1)
        self.walltimeminutes=np.int64(0)


    ################### 
    ## Create the job for QE pw.x calculation
    ###################
    def createjobQEpw(self,jobname,pwinput,pwoutput):
        # Inputs
        # jobname is the name of the job
        # pwinput is the path of the input file
        # pwoutput is the path of the output file

        print('')
        print(f'#########################')
        print(f'## Creating job {jobname}')

        # Open file
        f=open(jobname,'w')

        # Write heather of job
        f.write(f'#!/bin/bash\n\n')
        f.write(f'#PBS -N {self.name}\n')
        f.write(f'#PBS -l nodes={self.nodes}\n')
        f.write(f'#PBS -l ppn={self.ppn}\n')
        f.write(f'#PBS -q {self.queue}\n')
        f.write(f'#PBS -j oe\n')
        f.write(f'#PBS -l walltime={self.walltimehours:02d}:{self.walltimeminutes:02d}:00\n\n')
        
        f.write(f'cd $PBS_O_WORKDIR\n')
        f.write(f'module load espresso\n')
        f.write(f'mpirun -np $PBS_NP pw.x {pwinput} {pwoutput}\n')

        f.close()
        

