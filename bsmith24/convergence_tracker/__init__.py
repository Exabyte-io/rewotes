import sys
import os
import importlib

from convergence_tracker import settings
from convergence_tracker.utilities import general_utilities


class Convergence_Tracker:

    def __init__(self, program_name, *args, **kwargs):
        self.program_name = program_name
        self.kwargs = kwargs
        self.convergence_variables = self.kwargs['convergence_variables']
        self.jobs = [self.get_job(*args, **kwargs) for job in self.convergence_variables]
 

    def get_job(self, *args, **kwargs):
        """
        As of right now, this function returns an instance of the EspressoJob class.
        """
        if self.program_name in settings.PROGRAMS_REGISTRY:
            program_class = settings.PROGRAMS_REGISTRY[self.program_name]
        else:
            program_class = None 
        class_name = program_class.split('.')[-1]
        module_name = '.'.join(program_class.split('.')[0:-1])
        return getattr(importlib.import_module(module_name), class_name)(*args, **kwargs)


    def run_convergence_test(self):
        [job.run_job(index) for index, job in enumerate(self.jobs)]
        [job.update_job_state() for job in self.jobs]
        while any(state in [job.job_state for job in self.jobs] for state in ('Q', 'R')):
            [job.update_job_state() for job in self.jobs]
            print([job.job_state for job in self.jobs])
        [job.update_calculation_status(index) for index, job in enumerate(self.jobs)]
        if False in [job.is_finished for job in self.jobs]:
            print('ERROR: Cannot find the files pw.out')
            sys.exit(0)
        [job.update_convergence_property_value(index, 'total_energy') for index, job in enumerate(self.jobs)]
        total_energies = [job.convergence_property_value for job in self.jobs]
        convergence_results = general_utilities.is_converged(total_energies, tolerance=1.0)
        print('Quantum Espresso convergence calculation complete.')
        print('Convergence has been achieved:', convergence_results[0])
        print('Converged value of ecutwfc:', self.convergence_variables[convergence_results[2]])
