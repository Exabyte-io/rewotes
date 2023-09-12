import warnings
from MaterialPropertyPredictor import MPRLoader


loader = None

def test_loading():
    global loader
    warnings.filterwarnings("ignore", category=UserWarning)
    api_key_file = "api_key.txt"
    with open(api_key_file, "r") as f:
        api_key = f.read().strip()
        loader = MPRLoader()
        loader.load_data(api_key, chemsys=["Si", "Si-Ge", "Ge-Si", "Ge"])
    assert len(loader) > 0
    assert len(loader) == len(loader.get_model_inputs())
    assert len(loader) == len(loader.formulas)


def test_input_and_output_len():
    assert len(loader.get_model_inputs()) == len(loader.get_model_outputs())

def test_training_and_testing_splits():
    # test_size = 0.3
    train_x, train_y = loader.get_train_data()
    assert len(train_x) == len(train_y)

    test_x, test_y = loader.get_test_data()
    assert len(test_x) == len(test_y)

    training_length = len(train_x)
    testing_length = len(test_x)
    assert training_length > testing_length