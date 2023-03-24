from convergence_tracker import KpointConvergenceTracker
from ase.io import read


if __name__ == '__main__':
    path = 'data'
    atoms = read(f'{path}/POSCAR')
    tracker = KpointConvergenceTracker(atoms, 'etotal', 400)
    tracker.run_calcs()
    tracker.plot_results()


