# This script shows an example for the convergence of ecutwfc

import sys
sys.path.append('..')
from src.espresso_utilities import Espresso_Calculation
from src.general_utilities import General_Utilities

files = '../files'
ecutwfcs = ['18','20','22','24','26','28','30']
kpoints = ['8 8 8 0 0 0']

espresso_job = Espresso_Calculation(files, ecutwfcs, kpoints)
espresso_job.run_convergence_test()
