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
