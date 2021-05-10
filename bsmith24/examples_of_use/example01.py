import sys
sys.path.append('..')
from convergence_tracker import Convergence_Tracker

program_name = 'espresso'
kwargs = {
    'input_template':'../files/pw.in',
    'submit_template':'../files/job.pbs',
    'convergence_variables':[22,23,24,25]
    }

convergence_test = Convergence_Tracker(program_name, **kwargs)
convergence_test.run_convergence_test()

