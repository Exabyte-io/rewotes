from ase.calculators.vasp import Vasp 
from ase.io import read 
from convergence.convergence import KptsConvergenceTracker 
import argparse
import os

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--directory', help='Run directory with VASP POSCAR, KPOINTS, INCAR and POTCAR files', required=True)
    parser.add_argument(
        '-t', '--threshold', help='Threshold to meet for KPOINTS convergence; default is VASP total energy', type=float)
    parser.add_argument(
        '-c', '--calculation_type', help='Type of periodic structure calculation, i.e., VASP, Espresso, etc.', required=True)
    parser.add_argument(
        '-p', '--convergence_property', help='Property to converge', required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # My VASP submission details, change/alter for your specific setup 
    executable_path = '/home/rtrottie/programs/vasp/5.4.4.2019/bin/vasp_std'
    time = '1:00:00'
    partition = 'debug'
    nodes = 1
    account = 'custws'
    
    vasp_run_command = 'srun --time %s --partition %s --nodes %s --account %s %s' % (time, partition, str(nodes), account, executable_path)
    commands_dct = {'VASP': vasp_run_command}
    
    args = argument_parser()
    kpts_convergence = KptsConvergenceTracker(args.convergence_property, args.threshold)
    converged_kpoints, convergence_table = kpts_convergence.run(args.calculation_type, args.directory, run_command=commands_dct[args.calculation_type])
    print(convergence_table)
    print('Converged with %s at property difference threshold = %s' % (str(converged_kpoints), str(args.threshold)))
