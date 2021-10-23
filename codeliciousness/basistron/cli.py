# -*- coding: utf-8 -*-
import os
from argparse import ArgumentParser, Namespace

from basistron import model


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
        "--property",
        type=str,
        choices=set(
            model.CalculatedReferenceProperty.__members__.keys(),
        ).union(
            model.ExperimentalReferenceProperty.__members__.keys(),
        ),
        required=True,
        help="property against which basis set selection is evaluated",
    )
    parser.add_argument(
        "--value",
        type=float,
        help="reference property value (optional)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        help="allowed tolerance as % deviation from input property value",
    )
    parser.add_argument(
        "--regime",
        type=str,
        choices=list(model.ReferenceRegime.__members__.keys()),
        default=model.ReferenceRegime.experimental.value,
        help="use experimental or calculated data as benchmark"
    )
    return parser


def process_args(args: Namespace) -> model.Execution:
    """Perform initialization validation."""
    property = model.validate_property(args.regime, args.property)
    if not os.path.isfile(args.xyz_path):
        raise FileNotFoundError(args.xyz_path)
    with open(args.xyz_path, "r") as f:
        xyz_data = [
            ln.strip().split() for ln in f.readlines()[2:]
        ]
    return model.Execution(
        xyz_data=xyz_data,
        property=property,
        regime=args.regime,
        value=getattr(args, "value", None),
        tolerance=getattr(args, "tolerance", None),
    )


