import os

class Espresso_Calculation:
    """
    This class contains utilities / methods for making and performing calculations with the Quantum Espresso program.
    """

    def __init__(self, template_files_path, ecutwfcs, kpoints):
        self.template_files_path = template_files_path
        self.ecutwfcs = ecutwfcs
        self.kpoints = kpoints
        self.espresso_input_template = None
        self.espresso_submit_template = None


    # Make a map of the workflow ..
    def run_convergence_test(self):
        self.get_espresso_input_template()
        self.get_espresso_submit_template()
        #self.make_espresso_driver()
        #self.run_espresso_convergence_test()
        pass


    def get_espresso_input_template(self):
        """
        This function reads the espresso input template located in self.template_files_path
        and stores it as a list of strings via updating self.espresso_input_template.
        
        Returns:
            None, but updates self.espresso_input_template.
        """

        espresso_input_template_path = os.path.join(self.template_files_path, 'pw.in')
        with open(espresso_input_template_path) as template_file:
            self.espresso_input_template = template_file.readlines()


    def get_espresso_submit_template(self):
        """
        This function reads the espresso submit file template located in self.template_files_path
        and stores it as a list of strings via updating self.espresso_submit_template.
        
        Returns:
            None, but updates self.espresso_submit_template.
        """

        espresso_submit_template_path = os.path.join(self.template_files_path, 'job.pbs')
        with open(espresso_submit_template_path) as template_file:
            self.espresso_submit_template = template_file.readlines()
