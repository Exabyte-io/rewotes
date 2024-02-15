from pymatgen.io import pwscf
import subprocess
import os
'''Input and output handling
Formats supported: Quantum Espresso (PWscf input)
'''
#TODO: add additional format support

def load_job(path,job_type = "pwscf"):
    '''
    General function for loading a job from an input--handles various job types (e.g., 'pwscf', ...)
    :param path: str
    :param job_type: str
    '''
    job_type = job_type.lower() #to avoid capitalization input issues
    if (job_type=="pwscf"):
        job = QEJob(path)
    else:
        ValueError("Unsupported job_type "+job_type)
    return job
class Job():
    '''
    Base class for a file to run calculations
    '''
    def __init__(self,path):
        self.path = path
class QEJob(Job):
    '''
    Class for running Quantum Espresso calculations. Extends Job()
    '''
    def __init__(self,path):
        super().__init__(path)
        self.input = self.load(path) #PWInput object
        self.output = None
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
        '''Saves and runs self.input'''
        base, ext = os.path.splitext(self.path)
        if (output_path is None):output_path = base+".out"
        qe_process_str = prefix_str+" pw.x -i "+self.path+" > "+output_path
        #TODO: add parallelism
        subprocess.run(qe_process_str,shell=True)
        self.output = pwscf.PWOutput(output_path)


