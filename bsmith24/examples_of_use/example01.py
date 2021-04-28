# This script shows an example for the convergence of ecutwfc

import sys
sys.path.append('..')
from src.espresso_utilities import Espresso_Calculation

files = '../files'
ecutwfcs = ['28','30','32']
kpoints = ['4 4 4 0 0 0']

espresso_job = Espresso_Calculation(files, ecutwfcs, kpoints)
#espresso_job.perform_convergence_test()
