from ase.calculators.vasp import Vasp 
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.optimize import QuasiNewton 
import argparse
from copy import deepcopy
import sys
import os
import numpy as np

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--directory', help='Run directory with VASP POSCAR, KPOINTS, INCAR and POTCAR files', required=True)
    parser.add_argument(
        '-t', '--threshold', help='Threshold to meet for KPOINTS convergence; default is VASP total energy', type=float)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    ### My VASP submission details, change/alter for your specific setup ###
    executable_path = '/home/rtrottie/programs/vasp/5.4.4.2019/bin/vasp_std'
    time = '1:00:00'
    partition = 'debug'
    nodes = 1
    account = 'custws'
    run_command = 'srun --time %s --partition %s --nodes %s --account %s %s' % (time, partition, str(nodes), account, executable_path)
    ###

    args = argument_parser()
    poscar_path = os.path.join(args.directory, 'POSCAR')
    incar_path = os.path.join(args.directory, 'INCAR')
    kpoints_path = os.path.join(args.directory, 'KPOINTS')
    potcar_path = os.path.join(args.directory, 'POTCAR')

    aaa = AseAtomsAdaptor()
    pmg_atoms = Structure.from_file(poscar_path)
    ase_atoms = aaa.get_atoms(pmg_atoms)

    calc = Vasp(command=run_command, directory=args.directory)
    calc.read_incar(incar_path)
    calc.read_kpoints(kpoints_path)
    calc.read_potcar(potcar_path)

    all_kpoints = []
    all_energies = []
    num_kpts = 1

    while len(all_energies) < 2 or np.absolute(np.subtract(all_energies[-1], all_energies[-2])) > args.threshold:
        kpoints = [num_kpts, num_kpts, num_kpts] # Come up with a different way to get KPOINTS
        calc.set(kpts=kpoints)
        c_atoms = deepcopy(ase_atoms)
        c_atoms.set_calculator(calc)
        print('Running with KPOINTS = %s' % str(kpoints))
        try:
            optimizer = QuasiNewton(c_atoms, trajectory=os.path.join(args.directory, 'convergence.traj'))
            optimizer.run()
            all_kpoints.append(kpoints)
            all_energies.append(c_atoms.get_potential_energy()) # Make this modular, not just energy difference
        except IndexError:
            # Can occur if too few KPOINTS for tetrahedron method
            print('Bad calculation, check input files')
        num_kpts += 1

    final_ediff = np.around(np.absolute(np.subtract(all_energies[-1], all_energies[-2])), 5)
    final_kpoints = str(all_kpoints[-2])

    print('All kpoints: %s' % str(all_kpoints))
    print('All energies: %s' % str(all_energies))
    print('%s < %s threshold, converged with %s KPOINTs' % (str(final_ediff), str(args.threshold), str(final_kpoints)))
    
