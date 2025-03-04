"""Test the MaterialData class

This module contains tests for the MaterialData class in the data
loading module.

Example:
    To test the MaterialData class, run the following command:

    $ pytest tests/test_data_load.py
"""

import sys
from pathlib import Path
import pytest

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.data_load import MaterialData  # noqa: E402


@pytest.fixture(scope="module")
def api_key() -> str:
    """Read the API key from a file"""
    with open("api_key.txt", "r", encoding="utf-8") as f:
        api = f.read().strip()
    return api


@pytest.fixture(scope="module")
def data(api_key: str) -> MaterialData:
    """Return a MaterialData object"""
    return MaterialData(api_key, save=False, band_gap=(0.5, 0.55))


class TestMaterialData:
    def test_empty_init(self) -> None:
        """Expect a ValueError when no arguments are passed"""
        with pytest.raises(TypeError):
            MaterialData()

    def test_bad_init(self) -> None:
        """Expect a ValueError when an invalid argument is passed"""
        with pytest.raises(ValueError):
            MaterialData(24512)

    def test_init(self, data: MaterialData) -> None:
        """Expect a MaterialData object to be created"""
        assert data is not None

        # Check that the attributes are set correctly
        assert data.materials is None
        assert data.dataframe is None
        assert len(data) == 0

    def test_repr(self, data: MaterialData) -> None:
        """Expect the __repr__ method to return a string"""
        assert isinstance(repr(data), str)

    def test_get_materials(self, data: MaterialData) -> None:
        """Expect the material data to be fetched, but not saved"""
        materials = data.get_materials()

        # check that the materials are fetched
        assert materials is not None
        assert data.materials is not None
        assert len(materials) > 0

    def test_get_data(self, data: MaterialData) -> None:
        """Expect the material data to be fetched, cleaned, and not saved"""
        data._file_data = Path("temp/materials_data.hdf5")
        df = data.get_data()

        # check that the data is fetched and cleaned
        assert df is not None
        assert data.dataframe is not None
        assert len(data) > 0
        assert len(data) == len(df)
        assert isinstance(data.dataframe, pd.DataFrame)

        # check that the output file is not created
        assert not data._file_data.exists()
