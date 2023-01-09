from pkg import Downloader, TrainingDataManager, Model, Material

# create a Downloader object using your materialsproject.org API key
downloader = Downloader(my_api_key)

# download the materialsproject database to a local variable
materials = downloader.download()

# create an object to manage your training and testing sets
training_data_manager = TrainingDataManager(materials)

# create a neural network that is compatible with your training data
model = Model(training_data_manager)

# train the model using your training data
model.train(training_data_manager)

# create a material whose band gap you wish to predict
material = Material()
material.composition[8] = 1
material.composition[1] = 2
material.composition_reduced[8] = 1
material.composition_reduced[1] = 2

# predict the band gap for the new material
# the model generates multiple samples for the given input
# and returns the mean and sample standard deviation of the samples
band_gap_mean, band_gap_std = model.predict_new_material(material)

# band_gap_mean is a floating-point value, and it can be interpreted
# as the most like value for the band gap for the new material

# band_gap_std is a floating-point value, and it can be interpreted as
# the uncertainty in the prediction of the band gap
