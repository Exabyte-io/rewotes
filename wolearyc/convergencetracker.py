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
import json

from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints 

"""Class structure."""

def modify_pw_kpoints(pw_in_file, kpoints):
    """
    Modifies the kpoint specificiation in a pw.in.

    Args:
        pw_in (str): contents of a pw.in file
        kpoints (tuple): 3-tuple with k-point grid

    Returns:
        str
    """
    result = ''
    kpoints_flag = False
    for line in pw_in_file:
        if kpoints_flag == False:
            result += line
        else:
            result += f"{kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0"
            kpoint_flag = False
        if 'K_POINTS' in line:
            assert('automatic' in line)
            kpoints_flag = True

    return(result)


def gen_qe_workflow(input_file_path, kpoints, name):
    """
    Constructs and uploades a qe Mat3ra workflow using a pw.in file.

    Args:
        input_file_path (str): path to a pw.in input file (directories templated)
        kpoints (tuple): 3-tuple with k-point grid
        name (str): name of the workflow

    Returns:
        dict: the workflow
    """
    workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)

    # Generate new workflow for input file
    workflow_body = workflow_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
    workflow_body["name"] = name
    with open(input_file_path, "r") as f:
        content = modify_pw_kpoints(f, kpoints)
        workflow_body["subworkflows"][0]["units"][0]["input"][0]["content"] = content
        print(content)
    workflow = workflow_endpoints.create(workflow_body)
    return(workflow)


def gen_qe_job(input_file_path, kpoints):
    """ 
    Constructs and runs a qe job on Mat3ra.
    
    Args:
        input_file_path (str): path to a pw.in input file (directories templated)
        kpoints (tuple): 3-tuple with k-point grid

    Returns:
        (dict,dict): the workflow and job
    """
    job_endpoints = JobEndpoints(*ENDPOINT_ARGS)

    # A descriptive name for the job and workflow (hacky)
    name = f"{input_file_path.split('/')[-2]}_{kpoints}"

    workflow = gen_qe_workflow(input_file_path, kpoints, name)

    config = {
        "owner"   : {"_id": ACCOUNT_ID},
        "workflow": {"_id": workflow["_id"]},
        "name"    : name,
        "compute" : job_endpoints.get_compute(cluster = 'master-production-20160630-cluster-001.exabyte.io',
                                              ppn=1,queue='SR')
    }

    job = job_endpoints.create(config)
    return(workflow, job)

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

    def is_converged(self,job,ref_job):
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
workflow, job = gen_qe_job('Si/pw.in',(1,1,1))
print(job['_id'])
print(workflow['subworkflows'][0]['units'][0]['flowchartId'])

print(job)
                        