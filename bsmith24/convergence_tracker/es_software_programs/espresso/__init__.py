import os

from convergence_tracker.es_software_programs import Job


class EspressoJob(Job):

    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.convergence_variables = kwargs['convergence_variables']
        self.covergence_property = 'total energy'
        self.convergence_property_value = None
        self.job_output_file = 'pw.out'

    def edit_input_template(self, convergence_variable):
        for index, line in enumerate(self.input_template):            
            if 'ecutwfc' in line:
                self.input_template[index] = 'ecutwfc = '+str(convergence_variable)+'\n'

    def edit_submit_template(self):
        pass

    def update_convergence_property_value(self, index, property_to_get):
       espresso_output_file = os.path.join(str(index), self.job_output_file)
       with open(espresso_output_file) as output_file:
            filelines = output_file.readlines()
       if property_to_get == 'total_energy':
           for line in filelines:
               if '!    total energy' in line:
                   self.convergence_property_value = float(line.split()[4])*13605.662285    

    def update_calculation_status(self, index):
        espresso_output_file = os.path.join(str(index), self.job_output_file)
        with open(espresso_output_file) as output_file:
            filelines = output_file.readlines()
        for line in filelines:
            if 'JOB DONE' in line:
                self.is_finished = True

    def run_job(self, index):
        os.mkdir(str(index))
        os.chdir(str(index))
        self.edit_input_template(self.convergence_variables[index])
        with open('pw.in','w') as input_file:
            for line in self.input_template:
                input_file.write(line)
        with open('job.pbs', 'w') as submit_file:
            for line in self.submit_template:
                submit_file.write(line)
        self.submit_job()
        os.chdir('..')

