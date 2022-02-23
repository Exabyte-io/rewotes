#!/usr/bin/env python3

import argparse

from .bss import energy, get_homo_lumo_gap, optimize_basis_set


def main():
    """
    Run the CLI
    """
    parser = parser_setup()
    args = parser.parse_args()

    prop = args.property.lower()

    geom = read_geom(args.geom)

    if args.optimize:
        optimize(geom, prop, args)
    else:
        run_property(geom, prop, args)


def read_geom(geom: str = "geom.xyz") -> str:
    """
    Read an XYZ geometry
    """
    with open(geom) as f:
        next(f)
        next(f)
        return f.read()


def parser_setup():
    """
    Setup the argument parser
    """
    parser = argparse.ArgumentParser(description="Run a property computation on a geometry")
    parser.add_argument(
        "-g",
        "--geom",
        help="The file from which to get the geometry [%(default)s]",
        type=str,
        default="geom.xyz",
    )
    parser.add_argument(
        "-f",
        "--functional",
        help="The functional to use [%(default)s]",
        type=str,
        default="BP86",
    )
    parser.add_argument(
        "-b",
        "--basis_set",
        help="The basis set to use [%(default)s]",
        type=str,
        default="def2-svp",
    )
    parser.add_argument(
        "-p",
        "--property",
        help="Which property to return [%(default)s]",
        type=str,
        default="homo_lumo_gap",
    )
    parser.add_argument(
        "-t",
        "--target",
        help="Target value for selected property",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-o",
        "--optimize",
        help="Find the optimal basis_set [%(default)s]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--tolerance",
        help="Tolerance allowed when optimizing the property [%(default)s]",
        type=float,
        default=0.01,
    )

    return parser


def run_property(geom: str, prop: str, args):
    """
    Run property on the given geometry
    """
    target = args.target
    mf = energy(geom, functional=args.functional, basis_set=args.basis_set)
    if prop == "homo_lumo_gap":
        homo_lumo_gap = get_homo_lumo_gap(mf)
        error = homo_lumo_gap - target
        percent_error = error / homo_lumo_gap * 100
        print(
            f"""\
HOMO/LUMO Gap = {homo_lumo_gap:7.5f}
{target=:>7.5f}
{error=:>8.5f}
{percent_error=:>4.2f}%
"""
        )
    else:
        raise NotImplementedError(f"Property {prop} is not yet implemented")


def optimize(geom: str, prop: str, args):
    if prop == "homo_lumo_gap":
        property_function = get_homo_lumo_gap
    else:
        raise NotImplementedError(f"Property {prop} is not yet implemented")

    energy_kwargs = {
        "geom": geom,
        "functional": args.functional,
    }
    best = optimize_basis_set(
        energy_kwargs,
        property_function,
        target=args.target,
        tolerance=args.tolerance,
        stop_early=False,
    )
    for key, value in best.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
