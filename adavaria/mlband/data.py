import mlband.mp as MP
from .cgcnn.cgcnn.data import CIFData, DataLoader, collate_pool, get_train_val_test_loader
from pathlib import Path
import pandas as pd


elem_embedding_file='mlband/cgcnn/data/sample-regression/atom_init.json'

def get_list_of_materials(filter=True, num_chunks=1):
    fields=['material_id', 'composition', 'band_gap', 'structure']
    data = MP.mpr.materials.summary.search(
        # chemsys="Si-O", 
        fields=fields,
        chunk_size=1000,
        num_chunks=num_chunks,
    )

    df = pd.DataFrame({
        'material_id': [doc.material_id.__str__() for doc in data],
        'band_gap': [doc.band_gap for doc in data],
        'structure': [doc.structure for doc in data],
    })
    # df = df[fields]
    
    if filter:
        # Checking the missing values
        ind = ~(df['band_gap'] >= 0)
        print('Number of missing values: ', ind.sum())
        # Drop the missing values
        df = df.loc[~ind, :]
        # Reset the index
        df.reset_index(drop=True, inplace=True)
    
    return df, data

def create_dataset(df, path='data/cif_files'):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    for _, row in df.iterrows():
        structure = row['structure']
        material_id = row['material_id']
        # composition = row['composition']
        # band_gap = row['band_gap']
        filename = path / f'{material_id}.cif'
        writer = MP.CifWriter(structure)
        writer.write_file(filename)
        # print(f'Crystal structure {i+1} of {len(df)} written to {filename}.')
    
    df[['material_id', 'band_gap']].to_csv(path / 'id_prop.csv', index=False)
    import shutil
    shutil.copy(elem_embedding_file, path / 'atom_init.json')

