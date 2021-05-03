import subprocess

class Submit_Utilities:

    def __init__(self):
        pass


    @staticmethod
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


    @staticmethod
    def submit_driver_script():
        """
        A light-weight wrapper function to execute run.sh, which is the driver script
        for the calculations. This script makes the directory structure, and generates
        each espresso input file and pbs submission script.

        Returns:
            None, but executes the file run.sh.
        """

        subprocess.call(['sh', 'run.sh'])


    @staticmethod
    def get_pbs_job_states(job_ids):
        """
        This function gets the job status of pbs jobs via 'qstat' based on their job ids. It returns with 

        Args:
            job_ids (list of strings): A list containing the id for each job.

        Returns:
            job_states (list of strings): A list where each element is the state of the job

        """

        job_states = []
        queue_info = subprocess.check_output(['qstat']).decode('utf-8').split()
        for index, info in enumerate(queue_info):
            if info in job_ids:
                job_states.append(queue_info[index+5])    
        return job_states



class Job(Submit_Utilities):

    def __init__(self):
        super().__init__()
        self.job_id = None
        self.job_state = None
        self.is_submitted = None
        self.is_finished = None

