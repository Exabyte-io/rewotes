from .imports import *
import subprocess
import mlband.data
import torch
import torch.nn as nn
from . import utility


# cgcnn_path = Path('./cgcnn/').absolute().__str__()
# if cgcnn_path not in sys.path:
#     sys.path.append(cgcnn_path)


def train_by_original_code(config, ignore_warnings=True):
    """
    Training the model. Usage example:
    config = Config(epochs=50, learning_rate=0.005, batch_size=128)
    train(config)

    When using this function, cgccn must be cloned in the parent directory of mlband.
    git clone https://github.com/txie-93/cgcnn.git
    """
    print('Training the model...')
    results_dir = config.results_dir
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Add the cgcnn path to the system path
    cgcnn_path = Path('./cgcnn/').absolute().__str__()
    if cgcnn_path not in sys.path:
        sys.path.append(cgcnn_path)

    # ignore warnings
    if ignore_warnings:
        os.environ['PYTHONWARNINGS'] = "ignore"

    # Define the command and arguments
    # script_path = Path('./cgcnn/main.py').absolute().__str__()
    # script_path based on the relative path of this file
    script_path = Path(__file__).parent / 'cgcnn' / 'main.py'
    args = config.generate_args()
    command = ['python', script_path] + args

    # Start the process
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=results_dir)

    # Continuously showing the output
    while True:
        output = process.stdout.readline()
        if not output and process.poll() is not None:
            break
        if output:
            print(output.strip())

    # Close the stdout and waiting for the process to terminate
    process.stdout.close()
    return_code = process.wait()  # Capture the return code of the process

    # Check if the subprocess was successful
    if return_code != 0:
        raise Exception(f"Subprocess failed with return code {return_code}")

    print('Training and test evaluation completed!')


def train(config, ignore_warnings=True):
    """
    Training the model. Usage example:
    config = Config(epochs=50, learning_rate=0.005, batch_size=128)
    train(config)
    """
    print('Training the model...')
    args = config
    results_dir = config.results_dir
    Path(results_dir).mkdir(parents=True, exist_ok=True)

    if args.task == 'regression':
        best_mae_error = 1e10
    else:
        best_mae_error = 0.

    train_loader, val_loader, test_loader = mlband.data.get_data_loaders(config)
    normalizer = get_normalizer(config)
    model = get_cgcnn_model(config)
    criterion = nn.MSELoss()
    optimizer = get_optimizer(config, model)
    from torch.optim.lr_scheduler import MultiStepLR
    scheduler = MultiStepLR(optimizer, milestones=args.lr_milestones,
                            gamma=0.1)
    # Train the model
    for epoch in range(args.start_epoch, args.epochs):
        from .cgcnn.train import train_epoch, validate, save_checkpoint
        # train for one epoch
        train_epoch(config, train_loader, model, criterion, optimizer, epoch, normalizer)

        # evaluate on validation set
        mae_error = validate(config, val_loader, model, criterion, normalizer)

        if mae_error != mae_error:
            print('Exit due to NaN')
            sys.exit(1)

        scheduler.step()

        # remember the best mae_eror and save checkpoint
        if args.task == 'regression':
            is_best = mae_error < best_mae_error
            best_mae_error = min(mae_error, best_mae_error)
        else:
            is_best = mae_error > best_mae_error
            best_mae_error = max(mae_error, best_mae_error)
        save_checkpoint({
            'epoch': epoch + 1,
            'state_dict': model.state_dict(),
            'best_mae_error': best_mae_error,
            'optimizer': optimizer.state_dict(),
            'normalizer': normalizer.state_dict(),
            'args': vars(args)
        }, is_best, cd=results_dir)

    print('Training and test evaluation completed!')


def get_cgcnn_model(args,
                    orig_atom_fea_len=None,
                    nbr_fea_len=None,
                    atom_fea_len=None,
                    n_conv=None,
                    h_fea_len=None,
                    n_h=None,
                    ):
    from .cgcnn.data import CIFData
    from .cgcnn.model import CrystalGraphConvNet

    if orig_atom_fea_len is None or nbr_fea_len is None:
        dataset = CIFData(args.data_path)
        with utility.mute_warnings():
            structures, _, _ = dataset[0]
        orig_atom_fea_len = structures[0].shape[-1]
        nbr_fea_len = structures[1].shape[-1]

    model = CrystalGraphConvNet(orig_atom_fea_len, nbr_fea_len,
                                atom_fea_len=atom_fea_len or args.atom_fea_len,
                                n_conv=n_conv or args.n_conv,
                                h_fea_len=h_fea_len or args.h_fea_len,
                                n_h=n_h or args.n_h,
                                classification=False)
    if args.cuda:
        model.cuda()
    return model


def get_normalizer(config, ignore_warnings=True):
    from .cgcnn.data import CIFData, collate_pool
    from .cgcnn.utility import Normalizer
    from random import sample

    dataset = CIFData(config.data_path)

    with utility.mute_warnings(enabled=ignore_warnings):
        if len(dataset) < 500:
            warnings.warn('Dataset has less than 500 data points. '
                        'Lower accuracy is expected. ')
            sample_data_list = [dataset[i] for i in range(len(dataset))]
        else:
            sample_data_list = [dataset[i] for i in
                                sample(range(len(dataset)), 500)]
        _, sample_target, _ = collate_pool(sample_data_list)
        normalizer = Normalizer(sample_target)
    return normalizer


def load_model(config, orig_atom_fea_len=None, nbr_fea_len=None,):
    from .cgcnn.utility import Normalizer
    # Load the model
    best_checkpoint = torch.load(config.results_dir / 'model_best.pth.tar')
    model = get_cgcnn_model(config, orig_atom_fea_len=orig_atom_fea_len, nbr_fea_len=nbr_fea_len)
    model.load_state_dict(best_checkpoint['state_dict'])
    normalizer = Normalizer(torch.zeros(3))
    normalizer.load_state_dict(best_checkpoint['normalizer'])
    return model, normalizer


def evaluate_model(config):
    # Gettting the data loaders
    train_loader, val_loader, test_loader = mlband.data.get_data_loaders(config)

    # # Load the model
    # best_checkpoint = torch.load(config.results_dir/'model_best.pth.tar')
    # model = get_cgcnn_model(config)
    # model.load_state_dict(best_checkpoint['state_dict'])
    # normalizer = utility.Normalizer(torch.zeros(3))
    # normalizer.load_state_dict(best_checkpoint['normalizer'])
    # # criterion = nn.MSELoss()
    model, normalizer = load_model(config)

    # Evaluate the model
    for loader, name in zip([train_loader, val_loader, test_loader], ['train', 'val', 'test']):
        print('Evaluating the {} set...'.format(name))
        output_filename = os.path.join(config.results_dir, 'predictions_{}.csv'.format(name))
        df = predict_model(loader, model, normalizer, output_filename=output_filename)
        # print mse
        import sklearn.metrics
        mse = sklearn.metrics.mean_squared_error(df['True_Label'], df['Prediction'])
        mae = sklearn.metrics.mean_absolute_error(df['True_Label'], df['Prediction'])
        print(f'MSE for {name} set: {mse:.4f}')
        print(f'MAE for {name} set: {mae:.4f}')


def predict_model(loader, model, normalizer, output_filename=None, mute_warnings=True):
    model.eval()  # Set the model to evaluation mode
    predictions = []
    true_labels = []
    ids = []

    with torch.no_grad(), utility.mute_warnings(enabled=mute_warnings):
        for i, (inp, labels, batch_ids) in enumerate(loader):
            # Assuming your model and data are on the same device (CPU or CUDA)
            input_var = (inp[0], inp[1], inp[2], inp[3])

            # Compute output
            output = model(*input_var)
            denormed_output = normalizer.denorm(output.cpu())

            predictions.extend(denormed_output.view(-1).tolist())
            # Convert labels to scalar values
            if labels.dim() > 0:
                true_labels.extend(labels.squeeze().tolist())
            else:
                true_labels.append(labels.item())
            ids.extend(batch_ids)  # Assuming batch_ids is a list or similar
    
    df = pd.DataFrame({
        'ID': ids,
        'True_Label': true_labels,
        'Prediction': predictions
    })
    if output_filename is not None:
        # Save the predictions, true labels, and IDs to a CSV file
        print('Saving the predictions to {}'.format(output_filename))
        df.to_csv(output_filename, index=False)
    return df


def get_optimizer(config, model):
    import torch.optim as optim
    args = config
    if args.optim == 'SGD':
        optimizer = optim.SGD(model.parameters(), args.learning_rate,
                              momentum=args.momentum,
                              weight_decay=args.weight_decay)
    elif args.optim == 'Adam':
        optimizer = optim.Adam(model.parameters(), args.learning_rate,
                               weight_decay=args.weight_decay)
    else:
        raise NotImplementedError

    return optimizer


def predict(data_path, config=None, model=None, normalizer=None):
    if config is None:
        from .config import Config
        model_dir = Path(__file__).parent / 'pretrained_model'
        config = Config(data_path=data_path, results_dir=model_dir)

    data_loader = mlband.data.get_data_loader(config)
    if model is None or normalizer is None:
        model, normalizer = load_model(config)
    df = predict_model(data_loader, model, normalizer)
    return df
