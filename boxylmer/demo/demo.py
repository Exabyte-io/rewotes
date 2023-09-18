import sys
import os
demo_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(demo_dir) # feels very hacky just as a way to avoid pip install -e. 
sys.path.append(parent_dir)

from MaterialPropertyPredictor import MPRLoader, RandomForestBandGapModel, GradientBoostingBandGapModel
import matplotlib.pyplot as plt

from sklearn import metrics

def output_directory():
    script_dir = os.path.dirname(os.path.abspath(__file__))  
    output_dir = os.path.join(script_dir, "demo output") 

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir

def generate_parity_plot(y, pred, filename):

    rmse = metrics.mean_absolute_error(y, pred)  
    rmse = round(rmse, 2)
    plt.figure(figsize=(6, 6)) 
    plt.plot(y, pred, 'ro', markersize=8, markerfacecolor='red')
    max_value = max(max(y), max(pred))
    plt.plot([0, max_value], [0, max_value], 'k-', lw=2)

    plt.xlabel('y')
    plt.ylabel('pred')
    plt.xlim(0, max_value)
    plt.ylim(0, max_value)
    plt.gca().set_aspect('equal', adjustable='box') 
    plt.box(True)
    plt.grid(False) 
    plt.title(filename + " RMSE: " + str(rmse))
    
    plot_filename = os.path.join(output_directory(), filename + ".png")
    plt.savefig(plot_filename)

    assert os.path.isfile(plot_filename)

api_key_file = "api_key.txt"
with open(api_key_file, "r") as f:
    api_key = f.read().strip()

loader = MPRLoader(
    n_eigenvals=10
)
loader.load_data(
    api_key, 
    distance_method='fast',
    # elements=["Si", "Fe"],
    # chemsys=["Si-Ge", "Si", "Ge"],
    chemsys=["Fe", "O", "Si", "Fe-O", "Si-O", "Si-Fe", "Si-Fe-O"]
)

randomforest_model = RandomForestBandGapModel()
randomforest_model.fit(loader)
print(randomforest_model.model.feature_importances_)
y, pred = randomforest_model.parity(loader)
generate_parity_plot(y, pred, "random forest parity")

y, pred_train = randomforest_model.parity(loader, test_data_only=False)
generate_parity_plot(y, pred_train, "random forest parity - all data")



gradientboosting_model = GradientBoostingBandGapModel()
gradientboosting_model.fit(loader)
y, pred = gradientboosting_model.parity(loader)
generate_parity_plot(y, pred, "gradient boosting parity")

y, pred_train = gradientboosting_model.parity(loader, test_data_only=False)
generate_parity_plot(y, pred_train, "gradient boosting parity - all data")




# # Run this to identify better hyperparams
# randomforest_model.fit_hyperparameters(loader)
# y, pred = randomforest_model.parity(loader)
# generate_parity_plot(y, pred, "random forest parity - tuned hyperparameters")
# gradientboosting_model.fit_hyperparameters(loader)
# y, pred = gradientboosting_model.parity(loader)
# generate_parity_plot(y, pred, "gradient boosting parity - tuned hyperparameters")
