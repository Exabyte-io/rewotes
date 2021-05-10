import subprocess


def submit_pbs_job(job_script):
    """
    A light-weight wrapper function to the qsub tool using Python's 
    subprocess module.

    Args:
        job_script (string): The name of the 'pbs' file to be submitted
            using the call to 'qsub'.

    Returns:
        subprocess_output (string): Output from subprocess.check_output.
    """

    subprocess_output = subprocess.check_output(['qsub', job_script]).decode('utf-8').split()[0]
    return subprocess_output


def get_pbs_job_state(job_id):
    queue_info = subprocess.check_output(['qstat']).decode('utf-8').split()
    for index, info in enumerate(queue_info):
        if info in job_id:
            return queue_info[index+5]

