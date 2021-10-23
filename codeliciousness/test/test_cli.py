# -*- coding: utf-8 -*-

import pytest

from basistron import cli
from basistron.model import Execution

command = [
    "--property",
    "vibrational_frequency",
    "--value",
    "-100",
    "--xyz_path",
    "/path/to/file",
]


def mock_file(tmppath, h2):
    path = (tmppath / "h2.xyz").as_posix()
    with open(path, "w") as f:
        f.write(h2)
    return path


@pytest.mark.parametrize(
    "input_data, expected, raises",
    [
        ([], None, True),
        (
            ["--xyz_path", "/path/to/file"],
            None,
            True,
        ),
        (
            [
                "--xyz_path",
                "/path/to/file",
                "--property",
                "vibrational_frequency",
            ],
            {
                "xyz_path": "/path/to/file",
                "property": "vibrational_frequency",
            },
            False,
        ),
        (
            command,
            {
                "xyz_path": "/path/to/file",
                "property": "vibrational_frequency",
                "value": -100.0,
            },
            False,
        ),
    ],
)
def test_get_parser(input_data, expected, raises):
    parser = cli.get_parser()
    if raises:
        with pytest.raises(SystemExit):
            parser.parse_args(input_data)
    else:
        args = parser.parse_args(input_data)
        for key, val in expected.items():
            assert getattr(args, key) == val


def test_process_args(tmppath, h2, h2dat):
    path = mock_file(tmppath, h2)
    parser = cli.get_parser()
    args = parser.parse_args(command[:-1] + [path])
    ret = cli.process_args(args)
    assert isinstance(ret, Execution)
    assert ret.xyz_data == h2dat


def test_process_args_fail(tmppath, h2, h2dat):
    path = mock_file(tmppath, h2)
    parser = cli.get_parser()
    cmd = command.copy()
    args = parser.parse_args(cmd[:-1] + [path])
    args.property = "not_recognized"
    with pytest.raises(Exception):
        cli.process_args(args)
