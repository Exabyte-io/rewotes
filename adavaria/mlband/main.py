from mlband.imports import *
import mlband.data
import mlband.config
import mlband.model
import mlband.features


def main(num_chunks=1):
    data_set_path = 'data/mp-test-2/'
    results_dir = 'results/test-2/'

    print('Preparing the config...', flush=True)
    config = mlband.config.Config(
        data_path=data_set_path, 
        workers=4, 
        results_dir=results_dir,
        original_features=True,
        )

    print('Downloading the data from MP...' ,flush=True)
    df, data = mlband.data.get_list_of_materials(num_chunks=1)

    print('Creating the dataset...', flush=True)
    mlband.data.create_dataset(df=df, path=data_set_path)

    print('Creating the physical features...', flush=True)
    mlband.features.one_hot_encode_elements(config)
    
    print('Creating and training the model...', flush=True)
    mlband.model.train(config)

    print('Evaluating the training set...', flush=True)
    mlband.model.evaluate_model(config)


if __name__ == '__main__':
    main()
