from abc import ABC
from abc import abstractmethod

from convergence_tracker.utilities import job_utilities
from convergence_tracker.utilities import general_utilities


class Job(ABC):

    def __init__(self, *args, **kwargs):
        self.input_template = general_utilities.read_file(kwargs['input_template'])
        self.submit_template = general_utilities.read_file(kwargs['submit_template'])
        self.job_id = None
        self.is_submitted = False
        self.is_finished = False
        self.calculation_status = None
        self.convergence_property = kwargs['convergence_property']
        self.convergence_property_value = None
        self.output_file = None 

    @abstractmethod
    def edit_input_template(self):
        pass

    @abstractmethod
    def edit_submit_template(self):
        pass

    @abstractmethod
    def update_calculation_status(self):
        pass

    @abstractmethod
    def update_output_file_content(self, folder_name):
        pass

    @abstractmethod
    def update_convergence_property_value(self):
        pass

    @abstractmethod
    def run_job(self, folder_name):
        pass

    def submit_job(self) -> None:
        """
        A wrapper-like function to submit_pbs_job method in the job_utilities
        module. This function submits the job (to execute the electronic structure 
        program, obtains the job_id, and updates the is_submmited class variable.
        """
        self.job_id = job_utilities.submit_pbs_job('job.pbs')
        self.is_submitted = True    

    def update_job_state(self) -> None:
        """
        A wrapper-like function to get_pbs_job method in the job_utilities
        module. In this function, job_state is updated. 
        """
        self.job_state = job_utilities.get_pbs_job_state(self.job_id)

