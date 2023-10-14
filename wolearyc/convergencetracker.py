"""K-point convergence tracker."""

import sys
import json
import argparse
import time

import numpy as np

from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import get_jobs_statuses_by_ids
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints 

def wait_for_jobs_to_finish(endpoint: JobEndpoints, job_ids: list[str], poll_interval: int = 20):
    """
    Waits for jobs to finish and prints their statuses (adapted from Mat3ra API examples).
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.

    Args:
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to wait for
        poll_interval (int): poll interval for job information in seconds. Defaults to 10.
    """
    def spinning_cursor():
        while True:
            for cursor in '|/-\\':
                yield cursor
    cursor = spinning_cursor()
    
    counter = 0 
    statuses = get_jobs_statuses_by_ids(endpoint, job_ids)
    while True:
        if counter > poll_interval:
            counter = 0
            try:
                statuses = get_jobs_statuses_by_ids(endpoint, job_ids)
            except:
                # sometime the query fails. 
                continue
        errored_jobs   = len([status for status in statuses if status == "error"])
        active_jobs    = len([status for status in statuses if status == "active"])
        finished_jobs  = len([status for status in statuses if status == "finished"])
        submitted_jobs = len([status for status in statuses if status == "submitted"])
        print(f"{next(cursor)} [Mat3ra Jobs] Submitted:{submitted_jobs} Active:{active_jobs} Finished:{finished_jobs} Errored:{errored_jobs} (updates every {poll_interval} s)",
               end = '\r')
        time.sleep(0.5)
        counter += 0.5
        
        if all([status not in ["pre-submission", "submitted", "active"] for status in statuses]):
            print()
            break 

def modify_pw_kpoints(pw_in_lines: list[str], kpoints: tuple[int,int,int]):
    """
    Modifies the kpoint specificiation in a pw.in.

    Args:
        pw_in_lines (list): contents of a pw.in file
        kpoints (tuple): 3-tuple with k-point grid

    Returns:
        list[str]
    """
    result = []
    kpoints_flag = False
    for line in pw_in_lines:
        if kpoints_flag == False:
            result.append(line)
        else:
            result.append(f"{kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0")
            kpoint_flag = False
        if 'K_POINTS' in line:
            assert('automatic' in line)
            kpoints_flag = True

    return(result)

def modify_pw_ecut(pw_in_lines: list[str], ecutwfc: float, ecutrho: float):
    """
    Modifies the kinetic energy cutoffs in a pw.in. Adds an ecutrho entry if
    not specified.

    Args:
        pw_in_lines (list): contents of a pw.in file
        ecutwtf (float): wavefunction kinetic energy cutoff in Ry
        ecutrho (float): charge density/potential kinetic energy cutoff in Ry

    Returns:
        list
    """
    ecutrho_specified = 'ecutrho' in ''.join(pw_in_lines)
    result = []
    for line in pw_in_lines:
        if 'ecutwfc' in line:
            result.append(f"    ecutwfc = {ecutwfc}\n")
            if not ecutrho_specified:
                result.append(f"    ecutrho = {ecutrho}\n")
        elif 'ecutrho' in line:
            result.append(f"    ecutrho = {ecutrho}\n")
        else:
            result.append(line)

    return(result)

def gen_qe_workflow(input_file_path: str, kpoints: tuple[int,int,int], name: str,
                    ecutwfc: float, ecutrho: float):
    """
    Constructs and uploades a qe Mat3ra workflow using a pw.in file.

    Args:
        input_file_path (str): path to a pw.in input file (directories templated)
        kpoints (tuple): 3-tuple with k-point grid
        name (str): name of the workflow
        ecutwtf (float): wavefunction kinetic energy cutoff in Ry
        ecutrho (float): charge density/potential kinetic energy cutoff in Ry

    Returns:
        dict: the workflow
    """
    workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)

    # Generate new workflow for input file
    workflow_body = workflow_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
    workflow_body["name"] = name
    with open(input_file_path, "r") as f:
        content = modify_pw_kpoints(f.readlines(), kpoints)
        if ecutwfc and ecutrho:
            content = modify_pw_ecut(content, ecutwfc, ecutrho)
        content = ''.join(content)
        workflow_body["subworkflows"][0]["units"][0]["input"][0]["content"] = content
    workflow = workflow_endpoints.create(workflow_body)
    return(workflow)


def gen_qe_job(input_file_path: str, kpoints: tuple[int,int,int], material: dict, cores: int,
               ecutwfc: float, ecutrho: float):
    """ 
    Constructs and runs a qe job on Mat3ra.
    
    Args:
        input_file_path (str): path to a pw.in input file (directories templated)
        kpoints (tuple): 3-tuple with k-point grid
        material (dict): Mat3ra material
        cores (int): number of CPU cores to use
        ecutwtf (float): wavefunction kinetic energy cutoff in Ry
        ecutrho (float): charge density/potential kinetic energy cutoff in Ry

    Returns:
        (dict,dict): the workflow and job
    """
    job_endpoints = JobEndpoints(*ENDPOINT_ARGS)

    # A descriptive name for the job and workflow (hacky)
    name = f"{input_file_path.split('/')[-2]}_{kpoints}"

    workflow = gen_qe_workflow(input_file_path, kpoints, name, ecutwfc, ecutrho)

    config = {
        "owner"   : {"_id": ACCOUNT_ID},
        "workflow": {"_id": workflow["_id"]},
        "_material": {"_id": material["_id"]},
        "name"    : name,
        "compute" : job_endpoints.get_compute(cluster = 'master-production-20160630-cluster-001.exabyte.io',
                                              ppn=cores,queue='OR')
    }

    job = job_endpoints.create(config)
    return(workflow, job)

def find_next_kpoints(kpoints: tuple[int,int,int]):
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

def get_structure(pw_in_lines: list[str]):
    """
    Gets the structure from a pw_in string
    
    Args:
        pw_in_lines (list): contents of a pw.in file
    """
    elements = []; coordinates = []; cell_vectors = []

    # Hacky way of reading in atomic positions from pw.in
    # Unexpected spaces between ATOMIC_POSITIONS, CELL_PARAMETERS, and KPOINTS block
    # WILL lead to problems.
    reading_atomic_positions = False
    reading_cell_parameters = False
    for i,line in enumerate(pw_in_lines):        
        if 'ATOMIC_POSITIONS' in line:
            reading_atomic_positions = True
            continue
        elif 'CELL_PARAMETERS'in line:
            reading_atomic_positions = False
            reading_cell_parameters = True
            continue
        elif 'K_POINTS' in line:
            break

        if reading_atomic_positions:
            elements.append(line.split()[0])
            coordinates.append([float(n) for n in line.split()[1:]])
        elif reading_cell_parameters:
            cell_vectors.append([float(n) for n in line.split()])
        
    return(elements, coordinates, cell_vectors)

def get_cell_parameters(cell_vectors: list[list]):
    """Calculates cell parameters from a set of cell vectors.

    Args:
        cell_vectors(list): list of cell vectors

    Returns:
        float,float,float,float,float,float: a,b,c,alpha,beta,gamma
    """
    
    a = np.linalg.norm(cell_vectors[0])
    b = np.linalg.norm(cell_vectors[1])
    c = np.linalg.norm(cell_vectors[2])

    alpha = np.arccos(np.dot(cell_vectors[1], cell_vectors[2]) / (b * c))
    beta = np.arccos(np.dot(cell_vectors[0], cell_vectors[2]) / (a * c))
    gamma = np.arccos(np.dot(cell_vectors[0], cell_vectors[1]) / (a * b))
    alpha = np.degrees(alpha)
    beta = np.degrees(beta)
    gamma = np.degrees(gamma)

    return a, b, c, alpha, beta, gamma

def gen_material(pw_in_path: str):
    """
    Generates a Mat3ra material for the structure contained in a pw.in file

    Args:
        pw_in_path (str): path to pw.in file

    Returns:
        dict: the material
    """
    pw_in_lines = None
    with open(pw_in_path, 'r') as f:
        pw_in_lines = f.readlines()

    elements, coordinates, cell_vectors = get_structure(pw_in_lines)
    a, b, c, alpha, beta, gamma = get_cell_parameters(cell_vectors)
    
    config = {
        "name": pw_in_path.split('/')[-2],
        "basis": {
            "elements": [{"id":i+1, "value":element} for i,element in enumerate(elements)],
            "coordinates": [{"id":i+1, "value":vector} for i,vector in enumerate(coordinates)],
            "units": "crystal",
            "name": "basis",
            "cell" : cell_vectors,
        },
        "lattice": {
            "a": a,
            "b": b,
            "c": c,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "units": {"length": "angstrom", "angle": "degree"},
        },
        "tags": ["REST API"],
    }

    endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
    material = endpoint.create(config)
    print(f"    Mat3ra material id: {material['_id']}")
    return material

def print_job_info(kpoints: tuple[int,int,int], workflow: dict, job: dict):
    """
    Prints some info about a job.
    """
    print(f"           k-point grid: {kpoints[0]}x{kpoints[1]}x{kpoints[2]}")
    print(f"          Mat3ra job id: {job['_id']}")
    print(f"    Mat3ra flowchart id: {workflow['subworkflows'][0]['units'][0]['flowchartId']}")

def print_results(kpoints: tuple[int,int,int], energy: float, 
                  ref_kpoints: tuple[int,int,int], ref_energy: float):
    """
    Prints some info about a job.
    """
    print(f'Jobs complete.')
    print(f"    {energy} eV - {kpoints[0]}x{kpoints[1]}x{kpoints[2]}  VS")
    print(f"    {ref_energy} eV - {ref_kpoints[0]}x{ref_kpoints[1]}x{ref_kpoints[2]}")
    print(f"    Energy difference: {np.abs(energy-ref_energy)} eV")

class KConverger:
    """
    Superclass for k-point convergers.
    
    Args:
        input_file_dir (str): path to directory containing a pw.in file
        initial_kpoints (tuple): 3-tuple with initial k-point grid
        threshold (float): some positive threshold for convergence
        ecutwfc (float): wavefunction kinetic energy cutoff in Ry
        ecutrho (float): charge density/potential kinetic energy cutoff in Ry
        kpoints (list): all kpoints evaluated so far
        workflows (list): all workflows generated so far
        jobs (list): all jobs generated so far
        material (dict): the relevant material generated for all jobs 

    """

    def __init__(self,input_file_dir: str, threshold: float, initial_kpoints: tuple[int,int,int],
                 ecutwfc: float, ecutrho: float):
        self.input_file_dir = input_file_dir
        self.initial_kpoints = initial_kpoints
        self.threshold = threshold
        self.ecutwfc = ecutwfc
        self.ecutrho = ecutrho
        self.kpoints = []
        self.workflows = []
        self.jobs = []
        self.material = None

    def execute(self, cores: int):
        """ 
        Executes the convergence test.

        Args:
            cores(int): number of cores to use in calculations
        """
        all_kpoints = []
        all_workflows = []
        all_jobs = []
        job_endpoints = JobEndpoints(*ENDPOINT_ARGS)

        # Generate a material from the pw.in file
        print(f"Generating and uploading material...")
        self.material = gen_material(f"{self.input_file_dir}/pw.in")
        print(f"-" * 80)

        print(f"Running initial calculations...")
        # Generate initial two datapoints
        for kpoints in [self.initial_kpoints, find_next_kpoints(self.initial_kpoints)]:
            print(f'Generating and submitting job...')
            workflow, job = gen_qe_job(f'{self.input_file_dir}/pw.in', kpoints, self.material, cores,
                                       self.ecutwfc, self.ecutrho)
            job_endpoints.submit(job["_id"])
            print_job_info(kpoints, workflow, job)

            all_kpoints.append(kpoints)
            all_workflows.append(workflow)
            all_jobs.append(job)

        # Test for initial convergence
        
        wait_for_jobs_to_finish(job_endpoints, [j['_id'] for j in all_jobs], 10)
        converged, energy, ref_energy = self.check_convergence(all_workflows[-2], all_jobs[-2],
                                                               all_workflows[-1], all_jobs[-1])
        print_results(all_kpoints[-2], energy, all_kpoints[-1], ref_energy)
        if not converged:
            while not converged:
                print(f"-" * 80)
                kpoints = find_next_kpoints(kpoints)

                print(f'Generating and submitting job...')
                workflow, job = gen_qe_job(f'{self.input_file_dir}/pw.in', 
                                           kpoints, self.material, cores, 
                                           self.ecutwfc, self.ecutrho)

                job_endpoints.submit(job["_id"])
                print_job_info(kpoints, workflow, job)

                all_kpoints.append(kpoints)
                all_workflows.append(workflow)
                all_jobs.append(job)
                wait_for_jobs_to_finish(job_endpoints, [j['_id'] for j in all_jobs], 10)

                converged, energy, ref_energy = self.check_convergence(all_workflows[-2], all_jobs[-2],
                                                               all_workflows[-1], all_jobs[-1])
                print_results(all_kpoints[-2], energy, all_kpoints[-1], ref_energy)
                if converged:
                    print('Convergence reached.')

        # Store data in the object
        self.kpoints   = all_kpoints
        self.workflows = all_workflows
        self.jobs      = all_jobs

    def check_convergence(self,workflow: dict, job: dict, ref_workflow: dict, ref_job: dict):
        """ 
        Tests whether or not a job is converged.

        Args:
            workflow (dict): a workflow
            job (dict): job corresponding to workflow
            ref_workflow (dict): a more accurate workflow
            ref_job (dict): job corresponding to ref_workflow

        Returns:
            bool, float, float: whether converged, the score, and the reference score

        """
        raise NotImplementedError

class KEnergyConverger(KConverger):
    """
    Converges k-points with respect to total energy.
    """
    def __init__(self,input_file_dir: str, threshold: float, initial_kpoints: tuple[int,int,int],
                 ecutwfc: float, ecutrho: float):
        super().__init__(input_file_dir,threshold,initial_kpoints, ecutwfc, ecutrho)
    
    def check_convergence(self, workflow: dict, job: dict, ref_workflow: dict, ref_job: dict):
        """ 
        Checks for total convergence of a workflow/job pair

        Args:
            workflow (dict): a workflow
            job (dict): job corresponding to workflow
            ref_workflow (dict): a more accurate workflow
            ref_job (dict): job corresponding to ref_workflow

        Returns:
            bool, float, float: whether we have convergence, energy, and reference energy
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

def main(argv: list[str]):
    # Parse command line arguments
    parser = argparse.ArgumentParser(
                     prog='convergencetracker',
                     description='Finds a converged k-point grid for a Quantum Espresso calculation.',
                     epilog='Runs all QE calculations on Mat3ra.')
    parser.add_argument('input_file_dir', metavar='dir', type=str, 
                        help='path to directory containing a pw.in file')
    parser.add_argument('threshold', type=float,
                        help='convergence threshold (in eV)')
    parser.add_argument('-i', '--initialk', type=int, nargs = 3,
                        help='initial k-point grid')
    parser.add_argument('-e', '--ecuts', type=float, nargs = 2,
                        help='kinetic energy cutoffs (Ry) for wavefunctions and charge density/potential')
    parser.add_argument('-c', '--cores', type = int, 
                        help='number of cores to use in calculations')
    args = parser.parse_args(argv)
    
    # Default arguments:
    if args.initialk is None:
        args.initialk = (1,1,1)
    else:
        args.initialk = tuple(args.initialk)
    if args.cores is None:
        args.cores = 1
    if args.ecuts == None:
        args.ecuts = (None, None)
    
    # Print header
    print(f"=" * 80)
    print(f"CONVERGENCE TRACKER - by Willis O'Leary")
    print(f"")
    print(f"           Directory: {args.input_file_dir}")
    print(f"           Threshold: {args.threshold} eV")
    print(f"Initial k-point grid: {args.initialk[0]}x{args.initialk[1]}x{args.initialk[2]}")
    print(f"               Cores: {args.cores}")
    print(f"=" * 80)

    # Initialize and execute converger
    converger = KEnergyConverger(args.input_file_dir, args.threshold, 
                                 args.initialk, args.ecuts[0], args.ecuts[1])
    converger.execute(args.cores)

    result_kpoints = converger.kpoints[-2]

    # Print results
    print(f"=" * 80)
    print(f"CONVERGENCE REACHED")
    print(f"A {result_kpoints[0]}x{result_kpoints[1]}x{result_kpoints[2]} k-point grid is converged within {args.threshold} eV")
    print(f"    Required jobs: {len(converger.jobs)}")
    print(f"Finishing...")
    print(f"=" * 80)


if __name__ == "__main__":
   main(sys.argv[1:])
