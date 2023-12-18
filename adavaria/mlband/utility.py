import torch
import contextlib


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


@contextlib.contextmanager
def mute_warnings(enabled=True):
    """
    A context manager that temporarily suppresses warnings.
    Usage:
        with mute_warnings():
            warnings.warn('This will not be printed')
    """
    import warnings
    if enabled:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Add this block to configure logging for the pint library
            import logging
            pint_logger = logging.getLogger('pint')
            previous_level = pint_logger.getEffectiveLevel()
            pint_logger.setLevel(logging.ERROR)

            yield

            pint_logger.setLevel(previous_level)  # Reset the logging level
    else:
        yield

def mae(prediction, target):
    return torch.mean(torch.abs(target - prediction))

def mse(prediction, target):
    return torch.mean((target - prediction)**2)

def adjust_learning_rate(optimizer, epoch, k, lr):
    """Sets the learning rate to the initial LR decayed by 10 every k epochs"""
    assert type(k) is int
    lr = lr * (0.1 ** (epoch // k))
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr

def save_checkpoint(state, is_best, filename='checkpoint.pth.tar'):
    """Saves checkpoint to disk"""
    torch.save(state, filename)
    if is_best:
        import shutil
        shutil.copyfile(filename, 'model_best.pth.tar')

