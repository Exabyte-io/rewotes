import numpy as np


## THINGS TO FIX

class qepwinput:
    """This class describes the Quantum Espresso input for a pw.x calculation"""

    ################### 
    ## Initialization function for the class
    ###################
    def __init__(self):
        
        ## k-grid
        self.kgrid=np.array([1,1,1],dtype='int64')

    ################### 
    ## Load pw.x input file
    ## This function loads a pw.x input file into 
    ## the qepwinput class
    ###################
    def load(self,filename):
        # Inputs
        # filename is the path to pw.x input file

        print('')
        print(f'#########################')
        print(f'## Loading pw.x input from file {filename}')

        # Open file
        f=open(filename,'r')

        ## Retrieve k-grid
        # Loop over pw.x input file up to the kgrid information
        for line in f:
            x=line.split() 
            # if line is not empty
            if x:
                if (x[0]=='K_POINTS'): 
                    break
        # Read and save the k-grid information
        line=f.readline()
        x=line.split() 
        self.kgrid[0]=np.int64(x[0])
        self.kgrid[1]=np.int64(x[1])
        self.kgrid[2]=np.int64(x[2])

        f.close()
        

    ################### 
    ## Save pw.x input file
    ## This function saves an input file from a template
    ## replacing the current parameters of the qepwinput class
    ###################
    def save(self,template,filename):
        # Inputs
        # template is the path of the template pw.x input file
        # filename is the path to pw.x input file that will be created

        print('')
        print(f'#########################')
        print(f'## Saving pw.x input to file {filename}')

        # Open template file
        template=open(template,'r')

        # Open output file
        f=open(filename,'w')

        # Loop over pw.x input template and copy it up to 
        # the kgrid information
        for line in template:
            f.write(line)
            x=line.split() 
            # if line is not empty
            if x:
                if (x[0]=='K_POINTS'): 
                    break
        # write the k-grid information
        f.write(f'  {self.kgrid[0]} {self.kgrid[1]} {self.kgrid[2]} 0 0 0')

        f.close()
        template.close()






