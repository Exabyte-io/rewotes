from ase.calculators.vasp import Vasp 
from ase.io import read 
from convergence.tracker import ConvergenceHandler
import argparse
import os

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--directory', help='Run directory with VASP POSCAR, KPOINTS, INCAR and POTCAR files', required=True)
    parser.add_argument(
        '-t', '--threshold', help='Threshold to meet for KPOINTS convergence; default is VASP total energy', type=float)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # My VASP submission details, change/alter for your specific setup 
    executable_path = '/home/rtrottie/programs/vasp/5.4.4.2019/bin/vasp_std'
    time = '1:00:00'
    partition = 'debug'
    nodes = 1
    account = 'custws'
    run_command = 'srun --time %s --partition %s --nodes %s --account %s %s' % (time, partition, str(nodes), account, executable_path)

    args = argument_parser()
    ase_atoms = read(os.path.join(args.directory, 'POSCAR'))
    calc = Vasp(command=run_command, directory=args.directory)
    calc.read_incar(os.path.join(args.directory, 'INCAR'))
    calc.read_kpoints(os.path.join(args.directory, 'KPOINTS'))
    calc.read_potcar(os.path.join(args.directory, 'POTCAR'))
    starting_kpoints = [1, 1, 1]

    ch = ConvergenceHandler(ase_atoms, calc, starting_kpoints, args.threshold)
    converged_kpoints, convergence_table = ch.converge_kpoints()
    print(convergence_table)
    print('Converged with %s at property difference threshold = %s' % (str(converged_kpoints), str(args.threshold)))
