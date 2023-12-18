from .imports import *
import subprocess
import mlband.data
import torch
import torch.nn as nn
from . import utility


cgcnn_path = Path('./cgcnn/').absolute().__str__()
if cgcnn_path not in sys.path:
    sys.path.append(cgcnn_path)

def train(config, ignore_warnings=True):
    """
    Training the model. Usage example:
    config = Config(epochs=50, learning_rate=0.005, batch_size=128)
    train(config)
    """
    print('Training the model...')
    results_dir = config.results_dir
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Add the cgcnn path to the system path
    cgcnn_path = Path('./cgcnn/').absolute().__str__()
    if cgcnn_path not in sys.path:
        sys.path.append(cgcnn_path)

    #ignore warnings
    if ignore_warnings:
        os.environ['PYTHONWARNINGS'] = "ignore"

    # Define the command and arguments
    script_path = Path('./cgcnn/main.py').absolute().__str__()
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


def get_cgcnn_model(args, 
                    orig_atom_fea_len=None, 
                    nbr_fea_len=None,
                    atom_fea_len=None,
                    n_conv=None,
                    h_fea_len=None,
                    n_h=None,
                    ):
    from cgcnn.data import CIFData
    from cgcnn.model import CrystalGraphConvNet

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


def load_model(config):
    # Load the model
    best_checkpoint = torch.load(config.results_dir/'model_best.pth.tar')
    model = get_cgcnn_model(config)
    model.load_state_dict(best_checkpoint['state_dict'])
    normalizer = utility.Normalizer(torch.zeros(3))
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
        df = predict_model(loader, model, normalizer, output_filename)
        # print mse
        import sklearn.metrics
        mse = sklearn.metrics.mean_squared_error(df['True_Label'], df['Prediction'])
        print(f'MSE for {name} set: {mse:.4f}')


def predict_model(loader, model, normalizer, output_filename, mute_warnings=True):
    model.eval()  # Set the model to evaluation mode
    predictions = []
    true_labels = []
    ids = []

    with torch.no_grad(), utility.mute_warnings(enabled=mute_warnings):
        for i, (input, labels, batch_ids) in enumerate(loader):
            # Assuming your model and data are on the same device (CPU or CUDA)
            input_var = (input[0], input[1], input[2], input[3])

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

    # Save the predictions, true labels, and IDs to a CSV file
    print('Saving the predictions to {}'.format(output_filename))
    df = pd.DataFrame({
        'ID': ids,
        'True_Label': true_labels,
        'Prediction': predictions
    })
    df.to_csv(output_filename, index=False)
    return df


def get_optimizer(config, model):
    import torch.optim as optim
    args = config
    if args.optim == 'SGD':
        optimizer = optim.SGD(model.parameters(), args.lr,
                              momentum=args.momentum,
                              weight_decay=args.weight_decay)
    elif args.optim == 'Adam':
        optimizer = optim.Adam(model.parameters(), args.lr,
                               weight_decay=args.weight_decay)
        
    return optimizer

    

