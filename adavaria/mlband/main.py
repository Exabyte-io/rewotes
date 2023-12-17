from .imports import *
import mlband.data
import mlband.config
import mlband.model


def main(num_chunks=1):
    data_set_path = 'data/mp-test/'
    results_dir = 'results/test-1/'

    print('Downloading the data from MP...' ,flush=True)
    df, data = mlband.data.get_list_of_materials(num_chunks=1)

    print('Creating the dataset...', flush=True)
    mlband.data.create_dataset(df=df, path=data_set_path)

    print('Creating the config...', flush=True)
    config = mlband.config.Config(
        data_path=data_set_path, 
        workers=4, 
        results_dir=results_dir
        )
    
    print('Creating and training the model...', flush=True)
    mlband.model.train(config)

    print('Evaluating the training set...', flush=True)
