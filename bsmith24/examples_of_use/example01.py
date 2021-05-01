import sys
sys.path.append('..')
from src.espresso_utilities import Espresso_Calculation
from src.general_utilities import General_Utilities

files = '../files'
ecutwfcs = ['20','22','24','26','28']
kpoints = ['8 8 8 0 0 0']

espresso_job = Espresso_Calculation(files, ecutwfcs, kpoints)
#espresso_job.run_espresso_convergence_test()

print(espresso_job.total_energies)
espresso_job.get_each_espresso_total_energy()
print(espresso_job.total_energies)
convergence_results = General_Utilities.is_converged(espresso_job.total_energies, 1.0)
print(convergence_results)
