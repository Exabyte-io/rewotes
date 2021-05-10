from convergence_tracker.es_software_programs import Job

class EspressoJob(Job):


    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.covergence_property = 'total energy'


    def edit_input_template(self, convergence_variable):
        for index, line in enumerate(self.input_template):            
            if 'ecutwfc' in line:
                self.input_template[index] = 'ecutwfc = '+str(convergence_variable)+'\n'


    def edit_submit_template(self):
        pass


    def get_property(self, espresso_output_file, property_to_get):

       with open(espresso_output_file) as output_file:
            filelines = output_file.readlines()
       if property_to_get == 'total_energy':
           for line in filelines:
               if '!    total energy' in line:
                   return float(line.split()[4])*13605.662285


    def update_calculation_status(self, espresso_output_file):
        with open(espresso_output_file) as output_file:
            filelines = output_file.readlines()
        for line in filelines:
            if 'JOB DONE' in line:
                self.is_finished = True

