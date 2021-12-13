'''
    Materials Project doesn't appear to include HSE calculations,
    so we use data from Timur's paper, which has both GGA+HSE data.

    Draft structure:

    - Class (usable in other contexts) for encoding lattices as
    numerical vectors. Include compression option, i.e. with PCA
    for now (may improve learning). By default includes crystallographic
    and stoichiometric info of the material. Option to add 
    additional properties / info from intermediary calculations.

    - Class for ML that includes options for different models, e.g.
    regression, neural network, ... . Can choose loss function
    and specify hyperparameters depending on model chosen. Include
    error visualization and different cross-validation options.
'''

import numpy as np
from sklearn.decomposition import PCA

class Encode:
    '''
        Input needs to permit crystallographic info, stoichiometry,
        and additional features
    '''

    def __init__(self, X: np.array):
        '''Initialize data attributes'''
        self.dimX = np.shape(X)
        self.pca = self.pcaCompressed(X)

    def pcaCompressed(X):
        analyze = PCA(n_components=0.95, svd_solver='full')
        analyze.fit(X)
        return analyze.transform(X)

class Model:
    '''
        Needs to return learning model that can take input
        from Encode class. Allow to choose between multiple
        learning techniques as discussed.
    '''

    def __init__(self, model: str):
        '''
            model: string specifying learning technique from
                   fixed set of choices       
        '''

    def NN():
        return

    def regression():
        return