import torch
import numpy as np
from sklearn import metrics



class Normalizer(object):
    """
    Normalize a Tensor and restore it later. 
    Taken from: https://github.com/txie-93/cgcnn
    """

    def __init__(self, tensor=torch.zeros(3)):
        """tensor is taken as a sample to calculate the mean and std"""
        self.mean = torch.mean(tensor)
        self.std = torch.std(tensor)

    def norm(self, tensor):
        """Normalize a tensor"""
        return (tensor - self.mean) / self.std

    def denorm(self, normed_tensor):
        """Denormalize a normalized tensor"""
        return normed_tensor * self.std + self.mean

    def state_dict(self):
        """Returns a dictionary containing a whole state of the module."""
        return {'mean': self.mean,
                'std': self.std}

    def load_state_dict(self, state_dict):
        """Loads a state dict."""
        self.mean = state_dict['mean']
        self.std = state_dict['std']
        return self


class AverageMeter(object):
    """
    Computes and stores the average and current value.
    Taken from: https://github.com/txie-93/cgcnn
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

def class_eval(prediction, target):
    prediction = np.exp(prediction.numpy())
    target = target.numpy()
    pred_label = np.argmax(prediction, axis=1)
    target_label = np.squeeze(target)
    if not target_label.shape:
        target_label = np.asarray([target_label])
    if prediction.shape[1] == 2:
        precision, recall, fscore, _ = metrics.precision_recall_fscore_support(
            target_label, pred_label, average='binary')
        auc_score = metrics.roc_auc_score(target_label, prediction[:, 1])
        accuracy = metrics.accuracy_score(target_label, pred_label)
    else:
        raise NotImplementedError
    return accuracy, precision, recall, fscore, auc_score

def mae(prediction, target):
    return torch.mean(torch.abs(target - prediction))

def mse(prediction, target):
    return torch.mean((target - prediction)**2)