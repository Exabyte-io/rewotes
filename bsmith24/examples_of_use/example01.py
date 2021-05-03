import sys
sys.path.append('..')
from convergence_tracker.espresso_utilities import Espresso_Calculation

files = '../files'
ecutwfcs = ['22','23','24','25','26','27','28']
kpoints = ['8 8 8 0 0 0']

espresso_job = Espresso_Calculation(files, ecutwfcs, kpoints)
espresso_job.run_espresso_convergence_test()
