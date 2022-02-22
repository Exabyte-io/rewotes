#!/usr/bin/env python3

import argparse

from .bss import energy, get_homo_lumo_gap


def main():
    parser = argparse.ArgumentParser(description="Run a property computation on a geometry")
    parser.add_argument(
        "-g",
        "--geom",
        help="The file from which to get the geometry.",
        type=str,
        default="geom.xyz",
    )
    parser.add_argument(
        "-f", "--functional", help="The functional to use", type=str, default="BP86"
    )
    parser.add_argument(
        "-b", "--basis_set", help="The basis set to use.", type=str, default="def2-svp"
    )
    parser.add_argument(
        "-p",
        "--property",
        help="Which property to return [%s(default)]",
        type=str,
        default="homo_lumo_gap",
    )
    parser.add_argument("-t", "--target", type=float, help="Target value for selected property")

    args = parser.parse_args()
    target = args.target
    prop = args.property.lower()

    with open(args.geom) as f:
        next(f)
        next(f)
        geom = f.read()

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
        raise NotImplementedError("")


if __name__ == "__main__":
    main()
