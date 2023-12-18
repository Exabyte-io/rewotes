from .imports import *
from . import magpie

pt = magpie.get_magpie_data()


def one_hot_encode_elements(config=None, discrete_columns=None, continuous_columns=None, n_bins=10, filename=None):
    """
    One-hot encodes the elements. To get the list of available features, use the get_available_features() function.

    Parameters
    ----------
    config : Config
        The configuration object.
    discrete_columns : list of str
        A list of column names that contain discrete values.
    continuous_columns : list of str
        A list of column names that contain continuous values.
    n_bins : int
        The number of bins to use when binning continuous values.
    filename : str
        The name of the file to save the one-hot encoding dictionary to.
    """
    df = pt
    if config is not None:
        if config.original_features:
            import shutil
            # shutil.copy('mlband/atom_init.json', config.data_path)
            json_file = Path(Path(__file__).parent, 'atom_init.json')
            shutil.copy(json_file, config.data_path)
            # read the json file
            import json
            with open(json_file) as f:
                data = json.load(f)
            return data, 'original_features'
        filename = Path(config.data_path) / 'atom_init.json'
        discrete_columns = config.discrete_columns
        continuous_columns = config.continuous_columns

    if discrete_columns is None:
        discrete_columns = []
    if continuous_columns is None:
        continuous_columns = []

    encoded_df = pd.DataFrame()

    # One-hot encode discrete variables
    for col in discrete_columns:
        encoded_df = pd.concat([encoded_df, pd.get_dummies(df[col], prefix=col)], axis=1)

    # Bin and one-hot encode continuous variables
    for col in continuous_columns:
        bin_labels = [f"{col}_bin_{i}" for i in range(n_bins)]
        df[col + '_binned'] = pd.cut(df[col], bins=n_bins, labels=bin_labels)
        encoded_df = pd.concat([encoded_df, pd.get_dummies(df[col + '_binned'])], axis=1)

    # Create a dictionary with atomic numbers as keys and one-hot encoded features as values
    element_encoding = {str(row['atomic number']): encoded_df.iloc[idx].tolist()
                        for idx, row in df.iterrows()}
    # Converting T/F to 1/0
    element_encoding = {key: [int(i) for i in value] for key, value in element_encoding.items()}

    # Create a dictionary for feature names
    feature_names = {str(i): name for i, name in enumerate(encoded_df.columns)}

    if filename is not None:
        import json
        # Save the one-hot encoding dictionary as JSON
        with open(filename, 'w') as json_file:
            json.dump(element_encoding, json_file)

    return element_encoding, feature_names


def get_available_features():
    """
    Returns a list of available features.
    """
    return list(pt.columns)[1:]


if __name__ == '__main__':
    print(get_available_features())
