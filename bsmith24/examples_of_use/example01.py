import sys
sys.path.append('..')
from convergence_tracker import Convergence_Tracker

program_name = 'espresso'
kwargs = {
    'input_template':'../files/pw.in',
    'submit_template':'../files/job.pbs',
    'convergence_variable':'ecutwfc',
    'convergence_variables':[22,23,24],
    'convergence_property':'total_energy',
    'convergence_threshold':1e-5
    }

convergence_test = Convergence_Tracker(program_name, **kwargs)
convergence_test.run_convergence_test()

