import os

from convergence_tracker.es_software_programs import Job
from convergence_tracker.utilities import general_utilities


class EspressoJob(Job):


    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.convergence_variables = kwargs['convergence_variables']
        self.input_filename = 'pw.in'
        self.output_filename = 'pw.out'
        self.submit_filename = 'job.pbs'
        self.output_file_content = None


    def edit_input_template(self, convergence_variable) -> None:
        """
        This function updates self.input_template with the convergence variable

        Args:
            convergence_variable (int, float, or string): A particular list element of 
                self.convergence_variable. 
        """

        if self.convergence_property == 'ecutwfc':
            for index, line in enumerate(self.input_template):            
                if 'ecutwfc' in line:
                    self.input_template[index] = 'ecutwfc = '+str(convergence_variable)+'\n'
            

    def edit_submit_template(self) -> None:
        """
        This function updates self.submit_template with electronic structure program 
        specific information.
        """

        pass


    def update_output_file_content(self, index) -> None:
        """
        This function updates self.output_file with the content of the output
        file for a certain electronic structure program.
        """

        espresso_output_file = os.path.join(str(index), self.output_filename)
        self.output_file_content = general_utilities.read_file(espresso_output_file) 


    def update_convergence_property_value(self) -> None:
       """
       This function updates self.convergence_property_value with the
       value for the convergence property taken from the output file of a
       certain electronic structure program.
       """

       if self.convergence_property == 'total_energy':
           for line in self.output_file_content:
               if '!    total energy' in line:
                   self.convergence_property_value = float(line.split()[4])    


    def update_calculation_status(self) -> None:
        """
        This function updates self.is_finished according to if a condition is
        found within self.output_file_content.
        """

        for line in self.output_file_content:  
            if 'JOB DONE' in line:
                self.is_finished = True


    def run_job(self, index) -> None:
        """
        This function makes the directory structure and submits the jobs.
        """
        os.mkdir(str(index))
        os.chdir(str(index))
        self.edit_input_template(self.convergence_variables[index])
        with open(self.input_filename,'w') as input_file:
            for line in self.input_template:
                input_file.write(line)
        with open(self.submit_filename, 'w') as submit_file:
            for line in self.submit_template:
                submit_file.write(line)
        self.submit_job()
        os.chdir('..')
