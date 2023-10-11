"""Implementation of convergence tester."""

#Should this really be a package? A python script may be simpler.


#  User-facing interface.
#  User should be able to run the script directly (python ___.py) giving two 
#  arguments:
#     1. Path to input data (a pw.in or POSCAR, INCAR, KPOINTS, POTCAR)
#     2. Optionally, a kinetic energy cutoff (this will override the value in
#         INCAR)
#  
#  The program should print status messages, detailing status of submitted jobs, 
#  verifying job success, giving updates on the convergence, and finally
#  returning the result.

import sys

from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints 

"""Class structure."""

def run_qe_job(input_file_path, kpoints):
    """ 
    Constructs and runs a qe job on Mat3ra.
    
    Args:
        input_file_path (str): path to a pw.in input file (directories templated)
        kpoints (tuple): 3-tuple with k-point grid

    Returns:
        dict: the job
    """

class KConverger:
    """
    Superclass for k-point convergers.
    
    Args:
        input_file_dir (str): path to directory containing a pw.in file
        initial_kpoints (tuple): 3-tuple with initial k-point grid
        threshold (float): some positive threshold for convergence

    """

    def __init__(self,input_file_dir, initial_kpoints, threshold):
        self.input_file_dir = input_file_dir
        self.initial_kpoints = initial_kpoints
        self.threshold = threshold

    def execute(self):
        """ 
        Executes the convergence test.
        """

    def is_converged(self,job,ref_job)
        """ 
        Tests whether or not a job is converged.

        Args:
            job (dict): Mat3ra job
            ref_job (dict): reference Mat3ra job (should be more accurate than job)

        Returns:
            bool
        """
        raise NotImplementedError

class KEnergyConverger:
    """
    Converges k-points with respect to total energy.
    """
    def __init__(self,input_file_dir,initial_kpoints,threshold):
        super().__init__(input_file_dir,initial_kpoints,threshold)
    
    def is_converged(self,calculation,ref_calculation):
        raise NotImplementedError

# Below, we define the code 


                        