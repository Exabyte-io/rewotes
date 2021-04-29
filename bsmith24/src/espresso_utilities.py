import os

from src.general_utilities import General_Utilities


class Espresso_Calculation:
    """
    This class contains utilities / methods for making and performing calculations with the Quantum Espresso program.
    """

    def __init__(self, template_file_path, ecutwfcs, kpoints):
        self.template_file_path = template_file_path
        self.ecutwfcs = ecutwfcs
        self.kpoints = kpoints
        self.espresso_job_template = None


    # Make a map of the workflow ..
    def run_convergence_test(self):
        """
        This is a wrapper function for the workflow of the convergence test. It allows the user to
        perform the convergence test as: obj.run_convergence_test().

        Returns:
            None
        """
 
        self.get_espresso_job_template()
        self.update_espresso_job_template()       
        self.write_espresso_driver_script()


    def get_espresso_job_template(self):
        """
        A wrapper-like function for the get_job_template() function in the general_utilities module.
        
        Returns:
            None, but updates self.espresso_template.
        """
        self.espresso_job_template = General_Utilities.get_job_template(self.template_file_path)


    def update_espresso_job_template(self):
        """
        This function updates updates self.espresso_job_template with the kpoints and ecutwfcs
        values given by the user. 

        Returns:
            None, but updates self.espresso_job_template
        """
        assert self.espresso_job_template is not None
        self.espresso_job_template[1] = 'kpoints='+"'"+' '.join(self.kpoints)+"'"+'\n'
        self.espresso_job_template[2] = 'for ecut in '+' '.join(self.ecutwfcs)+'\n'


    def write_espresso_driver_script(self):
        """
        This function makes a file 'run.sh' with the contents of self.espresso_job_template.

        Returns:
            None
        """
        with open('run.sh', 'w') as file:
            for line in self.espresso_job_template:
                file.write(line)

