import warnings
from MaterialPropertyPredictor import MPRLoader
from MaterialPropertyPredictor import RandomForestBandGapModel
from MaterialPropertyPredictor import GradientBoostingBandGapModel


import os
import matplotlib.pyplot as plt
import pytest

def output_directory():
    output_dir = "test_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def generate_parity_plot(y, pred, filename):
    plt.figure(figsize=(6, 6)) 
    plt.plot(y, pred, 'ro', markersize=8, markerfacecolor='red')
    max_value = max(max(y), max(pred))
    plt.plot([0, max_value], [0, max_value], 'k-', lw=2, label='Parity Line')

    plt.xlabel('y')
    plt.ylabel('pred')
    plt.xlim(0, max_value)
    plt.ylim(0, max_value)
    plt.gca().set_aspect('equal', adjustable='box') 
    plt.box(True)
    plt.grid(False) 
    plt.title(filename)
    
    plot_filename = os.path.join(output_directory(), filename + ".png")
    plt.savefig(plot_filename)

    assert os.path.isfile(plot_filename)


loader = None
def test_loading():
    global loader
    warnings.filterwarnings("ignore", category=UserWarning)
    api_key_file = "api_key.txt"
    with open(api_key_file, "r") as f:
        api_key = f.read().strip()
        loader = MPRLoader()
        loader.load_data(
            api_key, 
            elements=["O", "Si", "Ge"]
            # chemsys=["Si", "Si-Ge", "Ge-Si", "Ge"],
        )
    assert len(loader) > 0
    assert len(loader) == len(loader.get_model_inputs())
    assert len(loader) == len(loader.formulas)

def test_input_and_output_len():
    assert len(loader) > 1
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

randomforest_model = None
def test_random_forest_bandgap_model():
    global randomforest_model
    randomforest_model = RandomForestBandGapModel()
    randomforest_model.fit(loader)
    randomforest_model.predict(loader)
    y, pred = randomforest_model.parity(loader)
    assert len(y) == len(pred)

    generate_parity_plot(y, pred, "random forest parity")

gradientboosting_model = None
def test_gradient_boosting_bandgap_model():
    global gradientboosting_model
    gradientboosting_model = GradientBoostingBandGapModel()
    gradientboosting_model.fit(loader)
    gradientboosting_model.predict(loader)
    y, pred = gradientboosting_model.parity(loader)
    assert len(y) == len(pred)

    generate_parity_plot(y, pred, "gradient boosting parity")