"""
This module provides a class for loading and processing material data from the
Materials Project API.

The `MaterialData` class allows users to retrieve material data from the
Materials Project API, process and clean the data, and perform operations
such as splitting the data into train and test sets.

Example:
    ```python
    from pathlib import Path
    from src.data_load import MaterialData

    # load all materials with band gap between 0 and 1000 eV
    data = MaterialData(api_key, band_gap=(0.0, 1000.0))
    data.get_data()

    # split data into train and test sets for band gap prediction
    x_train, x_test, y_train, y_test, _ = data.split_data()
    ```

Classes:
    MaterialData: A class for loading and processing material data from the
    Materials Project API.
"""

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from mp_api.client import MPRester


class MaterialData:
    """A class for loading and processing material data from the Materials
    Project API.

    Extra fields from the Materials Project can be added with the `fields`
    parameter in the constructor.
    External data for each material can be added with the `add_data_columns`
    method.

    Args:
        api_key (str): The API key for accessing the material data.
        fields (list, optional): The list of fields to retrieve from the
        material data. Defaults to None.
        save (bool, optional): Whether to save the loaded data. Defaults to
        True.
        **kwargs: Additional keyword arguments to be passed to the material
        data API.

    Attributes:
        api_key (str): The API key for accessing the material data.
        fields (list): The list of fields to retrieve from the material data.
        save (bool): Whether to save the loaded data.
        kwargs (dict): Additional keyword arguments to be passed to the
        material data API.
        materials (list): The loaded material data.
        dataframe (pd.DataFrame): The processed material data.
        _dir_output (Path): The output directory for saving the data.
        _file_data (Path): The file path for saving the data.

    Methods:
        __init__: Initializes the MaterialData object.
        __repr__: Returns a string representation of the MaterialData object.
        __len__: Returns the number of rows in the material data.
        _fetch_materials: Fetches the material data from the API.
        get_materials: Returns the loaded material data.
        get_data: Returns the processed material data.
        split_data: Splits the material data into train and test sets.
        add_data_columns: Adds additional columns to the material data.
        _extract_data: Extracts and cleans the material data.
        _encode_data: Encodes the categorical columns in the material data.
    """

    def __init__(self, api_key: str, fields: list = None, save: bool = True, **kwargs):
        """
        Initialize the DataLoad object.

        Parameters:
        - api_key (str): The API key for accessing the data.
        - fields (list): The list of fields to retrieve from the data.
        Defaults to a predefined list of fields.
        - save (bool): Flag indicating whether to save the data.
        Defaults to True.
        - **kwargs: Additional keyword arguments.

        Raises:
        - ValueError: If the API key is not provided.
        """
        self.api_key: str = api_key
        self.fields: list[str] = fields or [
            "material_id",
            "composition_reduced",
            "symmetry",
            "structure",
            "band_gap",
        ]
        self.save: bool = save
        self.kwargs: dict = kwargs

        if not isinstance(self.api_key, str) or not self.api_key:
            raise ValueError("API key must be provided")

        self.materials: list = None
        self.dataframe: pd.DataFrame = None

        self._dir_output: Path = Path("./data")
        self._file_data: Path = self._dir_output / "materials_data.hdf5"

    def __repr__(self) -> str:
        """Return a string representation of the MaterialData object.

        Returns:
            str: A string representation of the MaterialData object.
        """
        return (
            f"MaterialData(api_key={self.api_key}, fields={self.fields}"
            + f", kwargs={self.kwargs})"
        )

    def __len__(self) -> int:
        """Return the number of rows in the material data frame.

        Returns:
            int: The number of rows in the material data frame.
        """
        return len(self.dataframe) if self.dataframe is not None else 0

    def _fetch_materials(self) -> None:
        """Retrieve the material data from the Materials Project API."""
        with MPRester(self.api_key) as mpr:
            self.materials = mpr.materials.summary.search(
                fields=self.fields, **self.kwargs
            )

    def get_materials(self) -> list:
        """Return the loaded Material Project API data.

        If the data has not been loaded, it will be fetched from the API.

        Returns:
            list: Material Project data for each material.
        """
        if self.materials is None:
            self._fetch_materials()
        return self.materials

    def get_data(self) -> pd.DataFrame:
        """Return the processed and cleaned material data.

        If the data has been cached, it will be loaded from the file.
        Otherwise, the data will be fetched from the API, cleaned, and
        saved to the file.

        Returns:
            pd.DataFrame: Material data
        """

        # load data if it exists
        if self.dataframe is None and self._file_data.exists():
            self.dataframe = pd.read_hdf(self._file_data, key="data")
        elif self.dataframe is None:
            self._extract_data()
            self._encode_data()

        if self.save:
            self._dir_output.mkdir(exist_ok=True, parents=True)
            self.dataframe.to_hdf(self._file_data, key="data", mode="w")

        return self.dataframe

    def split_data(
        self, target: str = "band_gap", test_size: float = 0.2, seed: int = 42
    ) -> tuple:
        """Split the material data into train and test sets.

        Parameters:
        - target (str): The target column for prediction. Defaults to
        "band_gap".
        - test_size (float): The proportion of the data to include in the test
        set. Defaults to 0.2.
        - seed (int): The random seed for splitting the data. Defaults to 42.

        Returns:
        - tuple: A tuple containing the train and test sets of the input
        features and the target variable, as well as the material IDs, and the
        scaler used to scale the data.
        """
        if self.dataframe is None:
            self.get_data()

        # extract ID for later use
        mpid = self.dataframe["id"]

        x = self.dataframe.drop(columns=[target, "id"]).to_numpy(dtype=np.float64)
        y = self.dataframe[target].to_numpy(dtype=np.float64).reshape(-1, 1)

        # drop rows with NaN entries in x or y
        mask_x = np.isnan(x).any(axis=1)
        mask_y = np.isnan(y).flatten()
        mask = mask_x | mask_y
        x = x[~mask]
        y = y[~mask]
        mpid = mpid[~mask]

        # test/train split
        x_train, x_test, y_train, y_test = train_test_split(
            x, y, test_size=test_size, random_state=seed, shuffle=True
        )

        # scale the training and testing data
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x_train)
        x_test = scaler.transform(x_test)

        # save the scaler
        if self.save:
            scaler_filename = self._dir_output / "scaler.save"
            joblib.dump(scaler, scaler_filename)

        return x_train, x_test, y_train, y_test, mpid, scaler

    def add_data_columns(self, data: dict) -> None:
        """Add additional columns to the material data.

        Note that the data will be added to the existing data frame and it is
        assumed that the data is already encoded if necessary.

        Parameters:
        - data (dict): A dictionary of additional columns to add to the data.
        """
        if self.dataframe is None:
            self._extract_data()
            self._encode_data()

        self.dataframe = self.dataframe.assign(**data)

        if self.save:
            self.dataframe.to_hdf(self._file_data, key="data", mode="w")

    def _extract_data(self) -> pd.DataFrame:
        """Extract and clean the material data from the API into a DataFrame
        for analysis.

        Returns:
            pd.DataFrame: The cleaned material data.
        """
        if self.materials is None:
            self._fetch_materials()

        cleaned_data = []
        for doc in self.materials:
            # extract subset of symmetry data
            keys = ["crystal_system", "symbol", "point_group"]
            d = doc.symmetry.dict()
            symmetry = dict((k, d[k]) for k in keys)
            symmetry["crystal_system"] = symmetry["crystal_system"].value

            # extract subset of structure data
            lattice = doc.structure.lattice
            structure = {
                "a": lattice.a,
                "b": lattice.b,
                "c": lattice.c,
                "alpha": lattice.alpha,
                "beta": lattice.beta,
                "gamma": lattice.gamma,
                "density": doc.structure.density,
            }

            # combine dicts
            data = {
                **{"id": doc.material_id.split("()")[0]},
                **doc.composition_reduced.as_dict(),
                **symmetry,
                **structure,
                **{"band_gap": doc.band_gap},
            }
            cleaned_data.append(data)

        # convert list of dicts to pandas, and fill missing values in elements
        self.dataframe = pd.DataFrame(cleaned_data).fillna(0)
        return self.dataframe

    def _encode_data(self) -> pd.DataFrame:
        """Encode the categorical columns in the material data.

        Returns:
            pd.DataFrame: The encoded material data.
        """
        if self.dataframe is None:
            self._extract_data()

        # one-hot encoding for categorical columns
        self.dataframe = pd.get_dummies(
            self.dataframe,
            columns=["crystal_system", "point_group", "symbol"],
            drop_first=True,
        )

        return self.dataframe
