import numpy
import torch

class Model(torch.nn.Module):
    def __init__(self, training_manager):
        super(Model, self).__init__()
        self.linear1 = torch.nn.Linear(training_manager.total_feature_width - 1, 100)
        self.relu1 = torch.nn.ReLU()
        self.linear2 = torch.nn.Linear(100, 100)
        self.relu2 = torch.nn.ReLU()
        self.linear3 = torch.nn.Linear(100,1)
    def forward(self, x):
        x = self.linear1(x)
        x = self.relu1(x)
        x = self.linear2(x)
        x = self.relu2(x)
        x = self.linear3(x)
        return x
    def train(self, training_manager, epochs = 5):
        criterion = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(self.parameters())
        for epoch in range(epochs):
            loss_for_epoch = 0.0
            for data in training_manager.training:
                x, y = data
                optimizer.zero_grad()
                predictions = self(x)
                loss = criterion(predictions, y)
                loss.backward()
                optimizer.step()
                loss_for_epoch += loss.item()
            print("loss", loss_for_epoch)
    def test_standard_deviation(self, training_manager):
        prediction_deltas_zero = []
        prediction_deltas_nonzero = []
        for data in training_manager.testing:
            x, y = data
            predictions = self(x)
            for i in range(len(y)):
                target = float(y[i][0])
                actual = float(predictions[i][0])
                delta = actual - target
                if target == 0.0:
                    prediction_deltas_zero.append(delta)
                else:
                    prediction_deltas_nonzero.append(delta)
        prediction_deltas_zero = numpy.array(prediction_deltas_zero, dtype=numpy.double)
        prediction_deltas_nonzero = numpy.array(prediction_deltas_nonzero, dtype=numpy.double)
        factor = 1.0 / training_manager.data_range_encoder_array.encoders[0].factor
        zero_mean = numpy.mean(prediction_deltas_zero) * factor
        zero_std = numpy.std(prediction_deltas_zero) * factor
        nonzero_mean = numpy.mean(prediction_deltas_nonzero) * factor
        nonzero_std = numpy.std(prediction_deltas_nonzero) * factor
        return zero_mean, zero_std, nonzero_mean, nonzero_std
