import sys
import os
import subprocess
import time

from convergence_tracker.general_utilities import General_Utilities
from convergence_tracker.job_utilities import Job, Submit_Utilities


class Espresso_Calculation:
    """
    This class contains utilities / methods for making and performing calculations with the Quantum Espresso program.
    """

    def __init__(self, template_file_path, ecutwfcs, kpoints):
        self.template_file_path = template_file_path
        self.ecutwfcs = ecutwfcs
        self.kpoints = kpoints
        self.total_energies = None
        self.espresso_job_template = None       
        self.espresso_jobs = [Job() for ecut in ecutwfcs]               


    def get_espresso_job_template(self):
        """
        A light wrapper-like function for the get_job_template() function in the 
        general_utilities module.
        
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


    def get_each_espresso_total_energy(self):
        """
        This function loops through the different espresso directories, and 
        extracts the total energy.

        Returns:
            None, but updates self.total_energies
        """

        self.total_energies = []
        for ecut in self.ecutwfcs:
            espresso_output_file = os.path.join(ecut, 'pw.out')
            self.total_energies.append(self.get_espresso_total_energy(espresso_output_file))


    def get_espresso_total_energy(self, espresso_output_file):
        """
        This function extracts the total energy a quantum espresso output file.

        Args:
            espresso_output_file (string): The the quantum espresso output file.

        Returns:
            total_energy (float): The total energy, in units of meV. 
        """

        command = 'grep "! *total energy" '+espresso_output_file
        total_energy_Ry = float(subprocess.getoutput([command]).split()[4])
        total_energy_meV = total_energy_Ry * 13605.662285137
        return total_energy_meV


    def submit_espresso_jobs(self):
        """
        Goes into each espresso job folder and submits each quantum espresso job.

        Returns:
            None, but will submit each quantum espresso job.
        """

        for job, ecut in enumerate(self.ecutwfcs):
            os.chdir(ecut)
            self.espresso_jobs[job].job_id = self.espresso_jobs[job].submit_pbs_job('job.pbs')
            self.espresso_jobs[job].is_submitted = True
            self.espresso_jobs[job].is_finished = False
            os.chdir('..')


    def update_each_espresso_calculation_status(self):
        """
        Updates the status of each quantum espresso calculation.

        Returns:
            None, but updates self.espresso_jobs[job].is_finished for each job.
        """

        for job, ecut in enumerate(self.ecutwfcs):
            espresso_output_file = os.path.join(ecut, 'pw.out')
            self.espresso_jobs[job].is_finished = self.update_espresso_calculation_status(espresso_output_file)


    def update_espresso_calculation_status(self, espresso_output_file):
        """
        Updates the status of a certain quantum espresso calculation. Checks two things:
        1) Checks if the file pw.out exists.
        2) Checks if the calculation is finished.

        Returns:
            None.
        """

        if os.path.isfile(espresso_output_file):
            with open(espresso_output_file) as output_file:
               for line in output_file.readlines():
                   if 'JOB DONE' in line:
                       return True 
        else:
            return False 


    def update_espresso_job_states(self):
        """
        Updates the states of each espresso Job().job_state.

        Returns:
            None.
        """

        job_ids = [job.job_id for job in self.espresso_jobs]
        job_states = Submit_Utilities.get_pbs_job_states(job_ids)
        for index, state in enumerate(job_states):
            self.espresso_jobs[index].job_state = state


    def run_espresso_convergence_test(self):
        """
        This is a wrapper function for the workflow of the convergence test. It allows the user to
        perform the convergence test as: obj.run_convergence_test().

        Returns:
            None.
        """

        self.get_espresso_job_template()
        self.update_espresso_job_template()
        General_Utilities.write_driver_script(self.espresso_job_template)
        Submit_Utilities.submit_driver_script()
        self.submit_espresso_jobs()
        self.update_espresso_job_states()
        while any(state in [job.job_state for job in self.espresso_jobs] for state in('Q', 'R')): 
            #print([job.job_state for job in self.espresso_jobs])
            self.update_espresso_job_states()
            time.sleep(10)
        print('Calculations have finished. Checking for Quantum Espresso output files')
        self.update_each_espresso_calculation_status()
        if False in [job.is_finished for job in self.espresso_jobs]:
            print('ERROR: Cannot find the files pw.out')
            sys.exit(0)
        self.get_each_espresso_total_energy()
        convergence_results = General_Utilities.is_converged(self.total_energies, tolerance=1.0)
        print('Quantum Espresso convergence calculation complete.')
        print('Convergence has been achieved:', convergence_results[0])
        print('Converged value of ecutwfc:', self.ecutwfcs[convergence_results[2]])

