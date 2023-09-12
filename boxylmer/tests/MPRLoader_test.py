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