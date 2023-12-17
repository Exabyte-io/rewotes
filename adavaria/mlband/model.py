from .imports import *
import subprocess
import mlband.data
import torch
from cgcnn.model import CrystalGraphConvNet
import torch.nn as nn


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
    script_path = './cgcnn/main.py'
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
    process.wait()

    print('Training and test evaluation completed!')

def evaluate(config):
    # Gettting the data loaders
    train_loader, val_loader, test_loader = mlband.data.get_data_loaders(config)

    # Load the model
    best_checkpoint = torch.load('model_best.pth.tar')
    model = get_cgcnn_model(config)
    model.load_state_dict(best_checkpoint['state_dict'])
    normalizer = best_checkpoint['normalizer']
    criterion = nn.MSELoss()
    
    # Evaluate the model
    for loader, name in zip([train_loader, val_loader, test_loader], ['train', 'val', 'test']):
        print('Evaluating the {} set...'.format(name))
        predict(loader, model, criterion, normalizer)


def get_cgcnn_model(args):
    from cgcnn.data import CIFData
    dataset = CIFData(args.data_path)
    structures, _, _ = dataset[0]
    orig_atom_fea_len = structures[0].shape[-1]
    nbr_fea_len = structures[1].shape[-1]
    model = CrystalGraphConvNet(orig_atom_fea_len, nbr_fea_len,
                                atom_fea_len=args.atom_fea_len,
                                n_conv=args.n_conv,
                                h_fea_len=args.h_fea_len,
                                n_h=args.n_h,
                                classification=True if args.task ==
                                                       'classification' else False)
    if args.cuda:
        model.cuda()
    return model


def predict(loader, model, normalizer, output_filename):
    from torch.autograd import Variable
    model.eval()  # Set the model to evaluation mode
    predictions = []

    with torch.no_grad():
        for i, (input, _, _) in enumerate(loader):
            # Assuming your model and data are on the same device (CPU or CUDA)
            # Adjust the input handling here if your model expects a different format
            input_var = (Variable(input[0]), Variable(input[1]), input[2], input[3])

            # Compute output
            output = model(*input_var)
            denormed_output = normalizer.denorm(output.cpu())

            predictions += denormed_output.view(-1).tolist()

    # Save the predictions to a CSV file
    print('Saving the predictions to {}'.format(output_filename))
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for pred in predictions:
            writer.writerow([pred])

    return predictions

    

