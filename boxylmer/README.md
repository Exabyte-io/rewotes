## Project results
I would say this project was moderately successful given the constraints of a 5 day time limit and the concurrent management of graduate obligations. I had a lot of fun working on this and learned a bit along the way too! 

Two tree based models, that is, random forest and gradient boosting, were chosen due to their ease-of-implementation in scikit-learn and their reputation as good initial-choice models. Additionally, these models facilitated the 5 day time constraint as they are tolerant to lack of feature regularization. A graph convolution approach would likely be far superior in performance, however, a number of caveats to these kinds of models would have put this project out of a 5 day timeline: 
- Regularization of features is much more important, requiring additional time to test and validate. 
- Training can be excessively slow, requiring either a proxy model or ample resources. 
- Predictions may require much larger datasets to generalize, and memorization is a larger concern that would have risked losing time to.
- There are *far* more hyperparameters to worry about. 
- My previous work regarding graph-based approaches were implemented in Julia, not Python, and there would have been about a day or so lost to finding the relevant libraries (which I later found were NetworkX, DGL, Pytorch Geometric, and some others anyway) and learning how to use them appropriately.

The resulting RMSE when using datasets with a smaller (3-9) number of atoms (which is a limitation of the featurization choices used in the data preparation) ranges from 0.25 to 0.3, which I'm happy with considering SOTA literature models report around 0.14 to 0.2 as [the current best (2018)](https://pubs.acs.org/doi/10.1021/acs.chemmater.8b00686).

see the [plan, notes, and progress](./plan,%20notes,%20and%20progress.md) for a more detailed log of what I struggled with, found success in, and what I would pursue given more time.




# Running this package

## Setup
No GPU-acceleration is used. Simply install requirements.txt

An API key is needed to use MPRester. The demo and tests will search for the key in a file called "api_key.txt" in the root git folder (i.e., alongside `./boxylmer`). Before running, create that file and put the key in there.

## Running tests
A very minimal testfile is included using Pytest that should indicate if the environment is set up correctly and basic functionality is working, but a complete set of unit tests is not yet implemented.  

Navigate to ./rewotes and run `./rewotes> pytest`. Output will be sent to `./test_output`


## Running the Demo
The demo demonstrates usage and hyperparameter tuning (which is commented out by default) of two proof-of-concept models. 
By default, the demo assumes you're running Python from the git folder root, where the api_key.txt is stored. To run the demo, navigate to the root (asummed to be "rewotes") and run `./rewotes> python .\boxylmer\demo\demo.py`. An importance matrix for the random forest model will be printed. 