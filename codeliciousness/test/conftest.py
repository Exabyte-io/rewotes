# -*- coding: utf-8 -*-

import pytest

import pathlib


@pytest.fixture
def tmppath(tmpdir):
    return pathlib.Path(tmpdir)


@pytest.fixture
def h2():
    return """2

H 0.0 0.0 0.0
H 0.0 0.0 0.7"""


@pytest.fixture
def h2dat():
    return [
        ("H", 0.0, 0.0, 0.0),
        ("H", 0.0, 0.0, 0.7),
    ]
