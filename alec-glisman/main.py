"""Main script that trains models using the XGBoostModels class.

The main function in this script is `main()`, which is responsible for
executing the script. It follows the steps mentioned above and does not
return any value.

To run this script, execute the `main()` function.

Example:
    python main.py

Note: Before running the script, make sure to provide the API key in a file
named "api_key.txt" located in the same directory as this script.
"""

from pathlib import Path

from src.data_load import MaterialData
from src.models import XGBoostModels, NeuralNetModels


def main() -> None:
    """
    Main function that executes the script.

    This function performs the following steps:
      1. Reads the API key from a file.
      2. Loads data using the MaterialData class.
      3. Splits the data into training and testing sets.
      4. Trains models using the XGBoostModels class.
      5. Trains models using a Neural Network.
      6. Prints a completion message.

    Returns:
        None
    """
    file_path = Path(__file__).resolve().parent
    seed = 42

    # API key is not included in the code for security reasons
    with open(file_path / "api_key.txt", "r", encoding="utf-8") as f:
        api_key = f.read().strip()

    # Load data
    data = MaterialData(api_key, band_gap=(0.0, 10.0))
    x_train, x_test, y_train, y_test, _, _ = data.split_data(seed=seed)

    # Train models
    xgb = XGBoostModels(x_train, y_train, x_test, y_test, save=True)
    xgb.train_models(seed=seed)
    xgb.evaluate_model()
    nn = NeuralNetModels(x_train, y_train, x_test, y_test, save=True)
    nn.train_models(seed=seed)
    nn.evaluate_model()

    # Notify user that the script has finished
    print("Script completed successfully.")


if __name__ == "__main__":
    main()
