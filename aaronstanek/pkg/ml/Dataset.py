import numpy
import torch
from typing import Tuple

class Dataset(torch.utils.data.Dataset):
    def __init__(self, data_list):
        data_np = numpy.array(data_list)
        self._length = len(data_np)
        self.y = torch.from_numpy(data_np[:,0:1])
        self.x = torch.from_numpy(data_np[:,1:])
    def __len__(self) -> int:
        return self._length
    def __getitem__(self, index: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.x[index], self.y[index]
