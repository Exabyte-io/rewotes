#!/usr/bin/env python
from convergence_tracker import KpointConvergenceTracker, PwCutoffConvergenceTracker
from ase.io import read
import argparse
from utils import *
import logging


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
    parser.add_argument("--min-ka", default=10)
    parser.add_argument("--max-ka", default=80)
    parser.add_argument("--num-ka", default=10)
    logger = logging.getLogger("ConvergenceTrackingDriver")
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    args = parser.parse_args()
    path = args.path
    try:
        package = PeriodicDftPackages[args.package]
    except KeyError:
        logger.error(messages("unsupported_package").format(args.package))
        graceful_exit()

    try:
        convergence_property = ConvergenceProperty[args.convergence_property]
    except KeyError:
        logger.error(messages("unsupported_property").format(args.convergence_property))
        graceful_exit()

    try:
        convergence_parameter = ConvergenceParameter[args.convergence_parameter]
    except KeyError:
        logger.error(
            messages("unsupported_parameter").format(args.convergence_parameter)
        )
        graceful_exit()

    # switch on the convergence parameter for which convergence tracker to build

    if convergence_parameter == ConvergenceParameter.kpoints:
        tracker = KpointConvergenceTracker(
            min_ka=args.min_ka,
            max_ka=args.max_ka,
            num_ka=args.num_ka,
            workdir=args.workdir,
            path=path,
            convergence_property=convergence_property,
            package=package,
            eps=float(args.epsilon),
            logger=logger,
        )

    elif convergence_parameter == ConvergenceParameter.encut:
        # not implemented just as an example
        try:
            kpoints = [int(k) for k in args.kpoints[1:-1].split(",")]
            assert len(kpoints) == 3
        except (ValueError, AssertionError):
            logger.error(messages("invalid_kpts").format(args.kpoints))
            graceful_exit()

        tracker = PwCutoffConvergenceTracker(
            path, convergence_property, package=package, kpoints=kpoints
        )

    tracker.read_input()
    tracker.setup_calcs()
    tracker.run_calcs()
    tracker.show_results()
