from .imports import *

# Data table: https://github.com/ramv2/magpie_python/tree/master/lookup-data


def prepare_data():
    from matminer.utils.data import MagpieData

    mag = pd.DataFrame(MagpieData().all_elemental_props).reset_index().rename(columns={'index': 'symbol'})
    mag['atomic number'] = np.arange(1, len(mag) + 1).astype(int)
    mag.set_index('atomic number', inplace=True)
    # mag['atomic number'] = np.arange(1, len(mag) + 1).astype(int)

    # makedirs(Path('Data', 'magpie'), exist_ok=True)
    mag.to_csv(Path(Path(__file__).parent, 'magpie_data.csv'), index=False)
    return mag


def get_periodic_table():
    from pymatgen.core.periodic_table import Element

    # Retrieve all elements from the periodic table provided by pymatgen
    elements = [Element.from_Z(i) for i in range(1, 119)]

    # Define rare earth elements (Sc, Y, and lanthanides), these element get 0 as their group number
    rare_earths = [21, 39] + list(range(57, 72))  # Atomic numbers for Sc, Y, and lanthanides

    # Extract the properties of each element into a DataFrame
    df = pd.DataFrame({
        'symbol': [e.symbol for e in elements],
        'atomic number': [e.Z for e in elements],
        'atomic mass': [e.atomic_mass for e in elements],
        # 'group': [e.group for e in elements],
        'group': [0 if e.Z in rare_earths else e.group for e in elements],
        'period': [e.row for e in elements]  # 'row' attribute corresponds to the period
    })
    return df


def get_magpie_data(fill_missing=True):
    filename = Path(Path(__file__).parent, 'magpie_data.csv')
    if not filename.exists():
        df = prepare_data()
    else:
        df = pd.read_csv(filename)

    if fill_missing:
        # Replace NaN with 0
        df = df.fillna(0)

    pt = get_periodic_table()

    # Merge the two DataFrames
    df = pt.merge(df, on=['symbol'], how='inner')
    return df


if __name__ == '__main__':
    data = prepare_data()
    get_magpie_data()
    pass
