#!/usr/bin/env python
from convergence_tracker import KpointConvergenceTracker, PwCutoffConvergenceTracker
from ase.io import read
import argparse
from utils import *


if __name__ == "__main__":
    # read string input and convert everything to objects
    parser = argparse.ArgumentParser()
    parser.add_argument("path")  # positional
    parser.add_argument("-p", "--package", default="vasp")
    parser.add_argument("-d", "--workdir", default=None)
    parser.add_argument("-k", "--kpoints", default="[3,3,3]")
    parser.add_argument("-e", "--epsilon", default=1e-5)
    parser.add_argument("--convergence_property", default="etotal")
    parser.add_argument("--convergence_parameter", default="kpoints")

    args = parser.parse_args()
    path = args.path
    try:
        package = PeriodicDftPackages[args.package]
    except KeyError:
        print(f"{args.package} is not a supported periodic DFT package")
        graceful_exit()

    try:
        convergence_property = ConvergenceProperty[args.convergence_property]
    except KeyError:
        print(f"{args.convergence_property} is not a supported convergence property")
        graceful_exit()

    try:
        convergence_parameter = ConvergenceParameter[args.convergence_parameter]
    except KeyError:
        print(f"{args.convergence_parameter} is not a supported convergence parameter")
        graceful_exit()

    # switch on the convergence parameter for which convergence tracker to build

    if convergence_parameter == ConvergenceParameter.kpoints:
        tracker = KpointConvergenceTracker(
            workdir=args.workdir,
            path=path,
            conv_property=convergence_property,
            package=package,
            eps=float(args.epsilon),
        )

    elif convergence_parameter == ConvergenceParameter.encut:
        # not implemented just as an example
        try:
            kpoints = [int(k) for k in args.kpoints[1:-1].split(",")]
            assert len(kpoints) == 3
        except (ValueError, AssertionError):
            print(f"invalid kpoints supplied: {args.kpoints}")
            graceful_exit()

        tracker = PwCutoffConvergenceTracker(
            path, convergence_property, package=package, kpoints=kpoints
        )

    tracker.read_input()
    tracker.setup_calcs()
    tracker.run_calcs()
    tracker.show_results()
