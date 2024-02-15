from pymatgen.io import pwscf
import subprocess
import os
'''Input and output handling
Formats supported: Quantum Espresso (PWscf input)
'''
#TODO: add additional format support

class Job():
    '''
    Base class for a file to run calculations
    '''
    def __init__(self,path,to_load = False):
        self.path = path
class QEJob(Job):
    '''
    Class for running Quantum Espresso calculations. Extends Job()
    '''
    def __init__(self,path,to_load = False):
        super().__init__(path,to_load)
        self.input = self.load(path) #PWInput object
    @staticmethod
    def load(path):
        '''Loads PWInput object from filepath'''
        return pwscf.PWInput.from_file(path)
    def update_k_points(self,k_points):
        '''Updates kpoints_grid for PWInput'''
        k_points_tuple = tuple(k_points) #force type
        self.input.__setattr__("kpoints_grid",k_points_tuple)
    def save(self,path):
        '''Saves a Quantum Espresso input file for self.input at path'''
        self.input.write_file(path)
    def run(self,output_path = None,prefix_str = ""):
        #TODO: check in Linux bash
        '''Saves and runs self.input'''
        base, ext = os.path.splitext(self.path)
        if (output_path is None):output_path = base+".out"
        qe_process_str = prefix_str+" pw.x -i "+self.path+" > "+output_path
        #TODO: add parallelism
        subprocess.run(qe_process_str,shell=True)

