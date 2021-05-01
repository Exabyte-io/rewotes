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

        subprocess_output = subprocess.check_output(['qsub', job_script]).decode('utf-8')
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


class Job(Submit_Utilities):

    def __init__(self):
        super().__init__()
        self.job_id = None
        self.is_submitted = None
        self.is_finished = None

