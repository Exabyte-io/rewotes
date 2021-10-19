# -*- coding: utf-8 -*-
import os

from argparse import ArgumentParser, Namespace

from .model import Driver, Property


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description="Run basis set selection"
    )
    parser.add_argument(
        "--xyz_path",
        type=str,
        required=True,
        help="path to an XYZ file available on the filesystem",
    )
    parser.add_argument(
        "--target_property",
        type=str,
        required=True,
        help="property against which basis set selection is evaluated",
    )
    parser.add_argument(
        "--reference_value",
        type=float,
        required=True,
        help="reference property value (assumes atomic units)",
    )
    return parser


def process_args(args: Namespace) -> Driver:
    if not os.path.isfile(args.xyz_path):
        raise FileNotFoundError
    # load xyz data
    with open(args.xyz_path, "r") as f:
        xyz_data = [
            ln.strip().split() for ln in f.readlines()[2:]
        ]
    # validate target property
    if not Property.is_valid_property(args.target_property):
        raise Exception("unrecognized property")
    # optional tolerance (default defined in Driver)
    tol = getattr(args, "reference_tolerance", None)
    return Driver(
        xyz_data=xyz_data,
        target_property=args.target_property,
        reference_value=args.reference_value,
        reference_tolerance=tol,
    )


