"""K-point convergence tracker."""

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
import argparse

import numpy as np

from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import wait_for_jobs_to_finish
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints 

# Hacky: get api keys 

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

def find_next_kpoints(kpoints):
    """
    Returns a larger kpoint grid with the same kpoint ratios.

    The algorithm is extremely simple - surely there's a more clever way to do this! 
    In addition, this function should really take into account the lengths of the
    reciprocal lattice vectors - this could be improved.

    Args:
        kpoints(tuple): 3-tuple of kpoint subdivisions

    Returns:
        tuple: 3-tuple of kpoints
    """
    
    trial_kpoints = np.array(kpoints, dtype='float64')
    for i in range(1,100): # hard limit to prevent superdense k-point grids
        min_value = np.min(trial_kpoints)
        trial_kpoints *= (min_value + 1)/min_value
        if np.allclose(trial_kpoints / np.rint(trial_kpoints), 1, atol=1e-4):
            break
        
    return tuple([int(i) for i in trial_kpoints])

def print_job_info(kpoints, workflow, job):
    """
    Prints some info about a job.
    """
    print(f'Running job...')
    print(f"           k-point grid: {kpoints[0]}x{kpoints[1]}x{kpoints[2]}")
    print(f"          Mat3ra job id: {job['_id']}")
    print(f"    Mat3ra flowchart id: {workflow['subworkflows'][0]['units'][0]['flowchartId']}")

class KConverger:
    """
    Superclass for k-point convergers.
    
    Args:
        input_file_dir (str): path to directory containing a pw.in file
        initial_kpoints (tuple): 3-tuple with initial k-point grid
        threshold (float): some positive threshold for convergence

    """

    def __init__(self,input_file_dir, threshold, initial_kpoints):
        self.input_file_dir = input_file_dir
        self.initial_kpoints = initial_kpoints
        self.threshold = threshold
        self.kpoints = []
        self.workflows = []
        self.jobs = []

    def execute(self):
        """ 
        Executes the convergence test.
        """
        all_kpoints = []
        all_workflows = []
        all_jobs = []
        job_endpoints = JobEndpoints(*ENDPOINT_ARGS)

        print(f"=" * 80)
        print(f"Running initial calculations...")
        # Generate initial two datapoints
        for kpoints in [self.initial_kpoints, find_next_kpoints(self.initial_kpoints)]:
            workflow, job = gen_qe_job(f'{self.input_file_dir}/pw.in', kpoints)
            job_endpoints.submit(job["_id"])
            print_job_info(kpoints, workflow, job)

            all_kpoints.append(kpoints)
            all_workflows.append(workflow)
            all_jobs.append(job)

        # Test for initial convergence
        
        wait_for_jobs_to_finish(job_endpoints, [j['_id'] for j in all_jobs], 20)
        converged, energy, ref_energy = self.check_convergence(all_workflows[-2], all_jobs[-2],
                                                               all_workflows[-1], all_jobs[-1])
        if converged:
            print('converged!', energy, 'vs', ref_energy)

        #if not :
        #    while np.max(kpoints[-1]) < 2: # simple sanity check
        #        #workflow, job = gen_qe_job('Si/pw.in', kpoints[0])
        #        kpoints.append(find_next_kpoints(kpoints[-1]))
        #        energies.

        # Store data in the object
        self.kpoints   = all_kpoints
        self.workflows = all_workflows
        self.jobs      = all_jobs

    def check_convergence(self,workflow,job,ref_workflow,ref_job):
        """ 
        Tests whether or not a job is converged.

        Args:
            job (dict): Mat3ra job
            ref_job (dict): reference Mat3ra job (should be more accurate than job)

        Returns:
            bool, float, float

        """
        raise NotImplementedError

class KEnergyConverger(KConverger):
    """
    Converges k-points with respect to total energy.
    """
    def __init__(self,input_file_dir,threshold,initial_kpoints):
        super().__init__(input_file_dir,threshold,initial_kpoints)
    
    def check_convergence(self,workflow,job,ref_workflow,ref_job):
        """ 
        Checks for total convergence.
        """
        raw_property_endpoints = RawPropertiesEndpoints(*ENDPOINT_ARGS)

        def get_energy(t_job, t_workflow):
            unit_flowchart_id = t_workflow['subworkflows'][0]['units'][0]['flowchartId']
            energy = raw_property_endpoints.get_property(t_job['_id'], 
                                                         unit_flowchart_id, 
                                                         'total_energy')['data']['value']
            return energy
        
        energy = get_energy(job, workflow)
        ref_energy = get_energy(ref_job,ref_workflow)
        converged = np.abs(ref_energy - energy) < self.threshold
        return converged, energy, ref_energy 

def main(argv):
    # Parse command line arguments
    parser = argparse.ArgumentParser(
                     prog='convergencetracker',
                     description='Finds a converged k-point grid for a Quantum Espresso calculation.',
                     epilog='Runs all QE calculations on Mat3ra.')
    parser.add_argument('input_file_dir', metavar='dir', type=str, 
                        help='path to directory containing a pw.in file')
    parser.add_argument('threshold', type=float,
                        help='convergence threshold (in eV)')
    parser.add_argument('-i', '--initialk', metavar = 'K', type=int, nargs = 3,
                        help='initial k-point grid')
    args = parser.parse_args(argv)
    
    # Default arguments:
    if args.initialk is None:
        args.initialk = [1,1,1]
    
    print(f"=" * 80)
    print(f"CONVERGENCE TRACKER - by Willis O'Leary")
    print(f"")
    print(f"           Directory: {args.input_file_dir}")
    print(f"           Threshold: {args.threshold} eV")
    print(f"Initial k-point grid: {args.initialk[0]}x{args.initialk[1]}x{args.initialk[2]}")
    print(f"=" * 80)

    # Initialize converger
    converger = KEnergyConverger(args.input_file_dir, args.threshold, 
                                 args.initialk)
    converger.execute()

if __name__ == "__main__":
   main(sys.argv[1:])
