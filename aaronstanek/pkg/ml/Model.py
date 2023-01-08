import numpy
import torch
from .TrainingDataManager import TrainingDataManager
from typing import Tuple

class Model(torch.nn.Module):
    def __init__(self, training_manager):
        super(Model, self).__init__()
        self.linear1 = torch.nn.Linear(training_manager.total_feature_width - 1, 100)
        self.drop1 = torch.nn.Dropout(p=0.2)
        self.relu1 = torch.nn.ReLU()
        self.linear2 = torch.nn.Linear(100, 100)
        self.drop2 = torch.nn.Dropout(p=0.2)
        self.relu2 = torch.nn.ReLU()
        self.linear3 = torch.nn.Linear(100, 30)
        self.drop3 = torch.nn.Dropout(p=0.2)
        self.relu3 = torch.nn.ReLU()
        self.linear4 = torch.nn.Linear(30, 10)
        self.drop4 = torch.nn.Dropout(p=0.2)
        self.relu4 = torch.nn.ReLU()
        self.linear5 = torch.nn.Linear(10,1)
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.linear1(x)
        x = self.drop1(x)
        x = self.relu1(x)
        x = self.linear2(x)
        x = self.drop2(x)
        x = self.relu2(x)
        x = self.linear3(x)
        x = self.drop3(x)
        x = self.relu3(x)
        x = self.linear4(x)
        x = self.drop4(x)
        x = self.relu4(x)
        x = self.linear5(x)
        return x
    def train(self, training_manager: TrainingDataManager, epochs: int = 5) -> None:
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
    def predict_with_error(self, x: torch.Tensor, iteration_count: int = None) -> Tuple[numpy.ndarray, numpy.ndarray]:
        if iteration_count is None:
            iteration_count = 10
        prediction_set = []
        for i in range(iteration_count):
            prediction_set.append(self(x).detach().numpy())
        prediction_set = numpy.array(prediction_set)
        return numpy.apply_along_axis(numpy.mean, 0, prediction_set), numpy.apply_along_axis(numpy.std, 0, prediction_set)
    def test_with_error(self, training_manager: TrainingDataManager, iteration_count: int = None) -> Tuple[float, float, float]:
        prediction_z_scores = []
        for data in training_manager.testing:
            x, y = data
            prediction_means, prediction_stds = self.predict_with_error(x, iteration_count=iteration_count)
            for i in range(len(y)):
                prediction_z_scores.append(abs(float((prediction_means[i][0] - y[i][0]) / (prediction_stds[i][0] + 10**-30))))
        prediction_z_scores.sort()
        return (
            prediction_z_scores[int(len(prediction_z_scores) * 0.25)],
            prediction_z_scores[int(len(prediction_z_scores) * 0.50)],
            prediction_z_scores[int(len(prediction_z_scores) * 0.75)],
            )
    def test_standard_deviation(self, training_manager: TrainingDataManager) -> Tuple[numpy.floating, numpy.floating, numpy.floating, numpy.floating]:
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
