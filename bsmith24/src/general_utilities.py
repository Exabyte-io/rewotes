import os

class General_Utilities:
    """
    This class contains methods that perform tasks not specific to any particular
    electronic strcutre program.
    """

    def __init__(self):
        pass


    @staticmethod
    def get_job_template(template_file_path):
        """
        This function reads the job template located in self.template_files_path
        and stores it as a list of strings.
        
        Args:
            template_file_path (string): The path to where the job template file is stored.
                Ex) '../template_files'

        Returns:
            template_filedata (list of strings): A list where each element is a line from 
                the job template file.
        """

        job_template_path = os.path.join(template_file_path, 'run.sh')
        with open(job_template_path) as template_file:
            template_filedata = template_file.readlines()
        return template_filedata


    @staticmethod
    def write_driver_script(job_template):
        """
        This function makes a file 'run.sh' with the contents of job_template.

        Returns:
            None, but makes the 'ready-to-go' run.sh script.
        """

        assert job_template is not None
        with open('run.sh', 'w') as file:
            for line in job_template:
                file.write(line)


    @staticmethod
    def is_converged(values, tolerance):
        """
        This function determines if convergence of a value in values w.r.t
        the tolerance parameter has been met.

        Args:
            values (list of floats): The values with which we are intrested to
               see whether or not convergence has been achieved.

            tolerance (float): The tolerance for total energy convergence.
                By default, the convergence tolerance is 1 meV.

        Returns:
            is_converged (list of Boolean, float, and int): Boolean is True if 
                convergence is reached, False if convergence is not reached. The 
                variable that is float is the converged value and the variable that 
                is int is its index in values. 
                Ex) [True, 8356.3296, 4] or
                    [False, None, None]
        """

        prev_value = 0 
        for count, curr_value in enumerate(values):
            if abs(curr_value - prev_value) <= tolerance:
                return [True, values[count-1], count-1]
            else:
                prev_value = curr_value
        if count == len(values)-1:
            return [False, None, None]


    @staticmethod
    def submit_job():
        """
        A light wrapper-like function to execute run.sh.
  
        Returns:
            None, but executes run.sh
        """
        os.system('sh run.sh')

