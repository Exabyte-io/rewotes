from abc import ABC
from abc import abstractmethod

from convergence_tracker.utilities import job_utilities
from convergence_tracker.utilities import general_utilities

class Job(ABC):

    def __init__(self, *args, **kwargs):
        self.input_template = general_utilities.read_file(kwargs['input_template'])
        self.submit_template = general_utilities.read_file(kwargs['submit_template'])
        self.job_id = None
        self.job_state = None
        self.is_submitted = False
        self.is_finished = False

    @abstractmethod
    def edit_input_template(self):
        pass

    @abstractmethod
    def edit_submit_template(self):
        pass

    @abstractmethod
    def get_property(self):
        pass

    @abstractmethod
    def update_calculation_status(self):
        pass

    def setup(self, folder_name):
        pass 

    def submit_job(self, filename):
        self.job_id = job_utilities.submit_pbs_job(filename)
        self.is_submitted = True
     
    def update_job_state(self):
        self.job_state = job_utilities.get_pbs_job_state(self.job_id)

