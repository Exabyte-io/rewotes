import torch

class ConfigOld:
    def __init__(self, data_options=None, task='regression', disable_cuda=False, workers=0, 
                 epochs=30, start_epoch=0, batch_size=64, learning_rate=0.01, 
                 lr_milestones=[100], momentum=0.9, weight_decay=0, print_freq=10, 
                 resume='', train_ratio=None, train_size=None, val_ratio=0.1, 
                 val_size=None, test_ratio=0.1, test_size=None, optim='SGD', 
                 atom_fea_len=64, h_fea_len=128, n_conv=3, n_h=1, orig_atom_fea_len=92,
                 nbr_fea_len=76, lr=0.01):
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
        print_freq: int (frequency of printing training status per batch)
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
        lr: float (learning rate)
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
        self.lr = lr

        self.cuda = not self.disable_cuda and torch.cuda.is_available()


class Config:
    def __init__(self, data_path, disable_cuda=True, workers=0, epochs=30, 
                 batch_size=128, learning_rate=0.01, lr_milestones=None, momentum=0.9, 
                 weight_decay=0, print_freq=10, resume='', train_ratio=None, 
                 train_size=None, val_ratio=0.1, val_size=None, test_ratio=0.1, 
                 test_size=None, optim='Adam', atom_fea_len=64, h_fea_len=128, 
                 n_conv=3, n_h=1, task='regression',
                 results_dir='./results',
                 ):
        """
        data_path: str (path to data), e.g. './data/mp/'
        disable_cuda: bool (disable cuda), default: True
        workers: int (number of workers for data loader), default: 0
        epochs: int (number of epochs to train), default: 30
        batch_size: int (batch size for training and prediction), default: 128
        learning_rate: float (initial learning rate), default: 0.01
        lr_milestones: list of int (epochs to decrease learning rate), default: [100]
        momentum: float (SGD momentum), default: 0.9
        weight_decay: float (weight decay), default: 0
        print_freq: int (frequency of printing training status per batch), default: 10
        resume: str (path to checkpoint)
        train_ratio: float (ratio of training set), default: None
        train_size: int (size of training set), default: None
        val_ratio: float (ratio of validation set), default: 0.1
        val_size: int (size of validation set), default: None
        test_ratio: float (ratio of test set), default: 0.1
        test_size: int (size of test set), default: None
        optim: str (optimizer, SGD or Adam), default: Adam
        atom_fea_len: int (length of atom feature vector), default: 64
        h_fea_len: int (length of hidden vector), default: 128
        n_conv: int (number of convolutional layers), default: 3
        n_h: int (number of hidden layers), default: 1
        task: str (regression or classification), default: regression

        results_dir: str (path to save results), default: './results'
        """
        if lr_milestones is None:
            lr_milestones = [100]

        self.data_path = data_path
        self.disable_cuda = disable_cuda
        self.cuda = not self.disable_cuda and torch.cuda.is_available()
        self.workers = workers
        self.epochs = epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.lr_milestones = lr_milestones
        self.momentum = momentum
        self.weight_decay = weight_decay
        self.print_freq = print_freq
        self.resume = resume
        self.train_ratio = train_ratio
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
        self.task = task
        self.results_dir = results_dir

    def generate_args(self):
        args = [self.data_path]
        if self.disable_cuda:
            args.append('--disable-cuda')
        args.extend(['--workers', str(self.workers),
                     '--epochs', str(self.epochs),
                     '--batch-size', str(self.batch_size),
                     '--lr', str(self.learning_rate),
                     '--lr-milestones'] + [str(milestone) for milestone in self.lr_milestones])
        args.extend(['--momentum', str(self.momentum),
                     '--weight-decay', str(self.weight_decay),
                     '--print-freq', str(self.print_freq),
                     '--resume', self.resume,
                     '--optim', self.optim,
                     '--atom-fea-len', str(self.atom_fea_len),
                     '--h-fea-len', str(self.h_fea_len),
                     '--n-conv', str(self.n_conv),
                     '--n-h', str(self.n_h),
                     '--task', self.task])

        if self.train_ratio is not None:
            args.extend(['--train-ratio', str(self.train_ratio)])
        elif self.train_size is not None:
            args.extend(['--train-size', str(self.train_size)])

        if self.val_ratio is not None:
            args.extend(['--val-ratio', str(self.val_ratio)])
        elif self.val_size is not None:
            args.extend(['--val-size', str(self.val_size)])

        if self.test_ratio is not None:
            args.extend(['--test-ratio', str(self.test_ratio)])
        elif self.test_size is not None:
            args.extend(['--test-size', str(self.test_size)])

        return args


if __name__ == '__main__':
    print('Testing Config...')
    args = Config('./data/mp/')
    print(args.generate_args())
    print('Done!')
