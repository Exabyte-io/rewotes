# Load in relevant libraries, and alias where appropriate
import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as transforms

from torch.utils.data import Dataset, DataLoader
class Data(Dataset):
    def __init__(self,X_train,Y_train):
        self.X=torch.from_numpy(X_train).float()
        self.Y=torch.from_numpy(Y_train).float()
        self.len=self.X.shape[0]
    def __getitem__(self,index):      
        return self.X[index], self.Y[index]
    def __len__(self):
        return self.len

from mlbands.neuralnets import LeNet3D, LeNet5


def reshapeX(array,channels=1):
    return array.reshape((array.shape[0],1,*array.shape[1:]))

def reshapeY(array):
    if len(array.shape)==1:
        return array.reshape(-1,1)
    else: 
        return array

def reshapeXY(data,channels=1):

    X,Y = data
    # return [reshapeX(X,channels=1),Y]
    return [reshapeX(X,channels=1),reshapeY(Y)]


# Device will determine whether to run the training on GPU or CPU.
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


class Machine:
    def __init__(self, batch_size=1, num_classes=1, learning_rate=0.001, num_epochs=10 ):
        '''Machine [ learning ] class for bulk calculations

        Parameters
        ----------
        batch_size  : int
            batch size
        num_classes : int
            number of output ("Y") classes
        learning_rate : float
            learning rate for optimizer
        num_epochs : int
            number of training epochs
        neuralnet: < neuralnets.neural_network Class>
            neural network model (default: LeNet3D)
        cost : <torch.nn Object>
            cost / loss function for training default: MSELoss
        input_channels : int
            number of input channels to neural network for (3-D structure, band_gap) (X,Y) combination, it is =1
        '''
        self.batch_size = batch_size 
        self.num_classes = num_classes 
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs

        self.neuralnet = LeNet3D
        self.cost = nn.MSELoss()            # the loss function
        # self.cost = nn.CrossEntropyLoss()
        self.input_channels = 1


    def learn(self, trainset, testset):
        '''neural network training and testing (validation)

        Parameters
        ----------
        trainset : < Group.X, Group.y >
            training set of x and y values
        testset : < Group.X, Group.y >
            testing set of x and y values
        '''
        # Train the model
        trainset = reshapeXY(trainset,channels=self.input_channels)
        train_dataset=Data(*trainset)

        train_loader = torch.utils.data.DataLoader(dataset = train_dataset,
                                                batch_size = self.batch_size,
                                                shuffle = True)

        model = self.neuralnet(self.num_classes).to(device)


        #Setting the optimizer with the model parameters and learning rate
        optimizer = torch.optim.Adam(model.parameters(), lr=self.learning_rate)

        #this is defined to print how many steps are remaining when training
        total_step = len(train_loader)


        # Training the model
        for epoch in range(self.num_epochs):
            for i, (images, labels) in enumerate(train_loader):  
                images = images.to(device)
                labels = labels.to(device)
                
                #Forward pass
                outputs = model(images)
                loss = self.cost(outputs, labels)
                    
                # Backward and optimize
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                        
                if (i+1) % 400 == 0:
                    print ('Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}' 
                                .format(epoch+1, self.num_epochs, i+1, total_step, loss.item()))


        # Test the model
        # In test phase, we don't need to compute gradients (for memory efficiency)

        testset = reshapeXY(testset,channels=self.input_channels)

        test_dataset=Data(*testset)

        test_loader = torch.utils.data.DataLoader(dataset = test_dataset,
                                        batch_size = self.batch_size,
                                        shuffle = True)


        with torch.no_grad():
            correct = 0
            total = 0
            for images, labels in test_loader:
                images = images.to(device)
                labels = labels.to(device)
                outputs = model(images)
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()

                print('predicted: {} ground-truth: {}'.format(predicted,labels ))

            print('Accuracy of the network on the test data: {} %'.format(100 * correct / total))

