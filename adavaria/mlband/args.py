import torch

class Config:
    def __init__(self, data_options=None, task='regression', disable_cuda=False, workers=0, 
                 epochs=30, start_epoch=0, batch_size=64, learning_rate=0.01, 
                 lr_milestones=[100], momentum=0.9, weight_decay=0, print_freq=10, 
                 resume='', train_ratio=None, train_size=None, val_ratio=0.1, 
                 val_size=None, test_ratio=0.1, test_size=None, optim='SGD', 
                 atom_fea_len=64, h_fea_len=128, n_conv=3, n_h=1, orig_atom_fea_len=92,
                 nbr_fea_len=76):
        """
        data_options: list of paths to data
        task: 'regression' or 'classification'
        disable_cuda: bool (disable cuda)
        workers: int (number of workers for data loader)
        epochs: int (number of epochs to train)
        start_epoch: int (epoch to start training. useful on resume)
        batch_size: int (batch size for training and prediction)
        learning_rate: float (initial learning rate)
        lr_milestones: list of int (epochs to decrease learning rate)
        momentum: float (SGD momentum)
        weight_decay: float (weight decay)
        print_freq: int (frequency of printing training status)
        resume: str (path to checkpoint)
        train_ratio: float (ratio of training set)
        train_size: int (size of training set)
        val_ratio: float (ratio of validation set)
        val_size: int (size of validation set)
        test_ratio: float (ratio of test set)
        test_size: int (size of test set)
        optim: str (optimizer)
        atom_fea_len: int (length of atom feature vector)
        h_fea_len: int (length of hidden vector)
        n_conv: int (number of convolutional layers)
        n_h: int (number of hidden layers)
        orig_atom_fea_len: int (length of original atom feature vector)
        nbr_fea_len: int (length of neighbor feature vector)
        """
        self.data_options = data_options if data_options is not None else []
        self.task = task
        self.disable_cuda = disable_cuda
        self.workers = workers
        self.epochs = epochs
        self.start_epoch = start_epoch
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.lr_milestones = lr_milestones
        self.momentum = momentum
        self.weight_decay = weight_decay
        self.print_freq = print_freq
        self.resume = resume
        self.train_ratio = train_ratio if train_ratio is not None else 1 - val_ratio - test_ratio
        self.train_size = train_size
        self.val_ratio = val_ratio
        self.val_size = val_size
        self.test_ratio = test_ratio
        self.test_size = test_size
        self.optim = optim
        self.atom_fea_len = atom_fea_len
        self.h_fea_len = h_fea_len
        self.n_conv = n_conv
        self.n_h = n_h
        self.orig_atom_fea_len = orig_atom_fea_len
        self.nbr_fea_len = nbr_fea_len

        self.cuda = not self.disable_cuda and torch.cuda.is_available()


args = Config()

