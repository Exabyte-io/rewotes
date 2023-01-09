import numpy
import torch
from typing import List, Tuple


class Dataset(torch.utils.data.Dataset):
    """Stores a collection of input rows and the associated target output for
    each row."""

    def __init__(self, data_list: List[numpy.ndarray]):
        """Create a Dataset instance from a list of 1D numpy arrays.

        Each numpy array is a row, with the first element being the
        target output for the row, and the remaining entries being the
        features.
        """
        data_np = numpy.array(data_list)
        self._length = len(data_np)
        self.y = torch.from_numpy(data_np[:, 0:1])
        self.x = torch.from_numpy(data_np[:, 1:])

    def __len__(self) -> int:
        """Return the number of rows in the Dataset."""
        return self._length

    def __getitem__(self, index: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """Return the nth row, as a tuple.

        The first element of the tuple is a torch Tensor of features.
        The second element is a singleton torch Tensor with the target
        value.
        """
        return self.x[index], self.y[index]
