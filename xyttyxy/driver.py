#!/usr/bin/env python
from convergence_tracker import KpointConvergenceTracker, PwCutoffConvergenceTracker
from ase.io import read
import argparse
from utils import *
import logging


if __name__ == "__main__":
    # read string input and convert everything to objects
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path to structure files (POSCAR-like, cif, extended-xyz) and calculation input files (only VASP INCAR-like files supported)")  # positional
    parser.add_argument("-p", "--package", default="vasp", help="periodic dft package to use. Only VASP is supported.")
    parser.add_argument("-d", "--workdir", default=None, help="where to put the output file (the convergence plot)")
    parser.add_argument("-k", "--kpoints", default="[3,3,3]", help="regular kpoints mesh for converging cutoff. not implemented yet. ")
    parser.add_argument("-e", "--epsilon", default=1e-5, help="convergence criteria, e.g. eV/Ry for VASP/QE, ")
    parser.add_argument("--convergence_property", default="etotal", help="the property to converge. Only total energy (etotal) is supported. ")
    parser.add_argument("--convergence_parameter", default="kpoints", help="calculation parameter to vary. Only kpoints is supported. ")
    parser.add_argument("--min-ka", default=10, help="for kpoints convergence. The minimum spacing, in <lattice vector length>*<number of kpoints>, assuming uniform spacing in all 3 dimensions")
    parser.add_argument("--max-ka", default=80, help="for kpoints convergence. The maximum spacing")
    parser.add_argument("--num-ka", default=10, help="for kpoints convergence. Number of k-point meshes to try. Note that actual number might be lower since identical meshes may be produced")
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
