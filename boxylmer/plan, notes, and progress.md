# Minimum Goals
- Create a model which will automatically fit and predict material band-gaps given an input. 

# The Plan
## Monday (primarily planning, familiarize with literature and MPRester)
### Day plan
- Create package structure
- Outline classes and necessary abstractions
    - Abstract class to load / curate data (DataLoader)
    - Abstract class to wrap prediction (BandGapPredictor)
- Fill out the data loader and curator
    - How would other kinds of data be incorporated?
        - Probably a dict that holds extra vectors and a validation function to check that every entry is the right corresponding length.
        - Depending on whether or not we can expect the inputs to be standardized in pymatgen's structure format, this may need to change a lot.
    - How would data from other sources be parsed? -> Another concrete DataLoader
    - Parsed data held in this way will be slightly inefficient, even for Python's standards, as np.arrays allocate.
    - But I think the result will be much more organized and it looks like structure.get_distance() will be the bottleneck anyhow.
    - Do we need to cache the data? -> Not within scope
        - Leave reading/writing method spaces for future dev if necessary
- Will leave relevant and useful literature sources as comments for reference later.

## Tuesday (Finish up data loader, begin model class)
### Day plan
- Identify appropriate structure for models and fill out the BandGapPredictor
    - in this limited case the input will be flat, but a 3D structure such as graphs will be more complicated and likely require a different representation or accessor method
    - The loader should hold onto the structure objects, so graph construction is still possible without weird methods / a class rewrite.
- Models: RandomForest (doable), GradientBoosting (doable), Graph Convolution (reach, likely won't have enough time this week)
- Padding -> Decide how to allocate it, how will the max size be determined and stored?
- Possibly, we could use an autoencoder to get over possible padding issues, if the tree based methods don't like the zero padding. 
    - Pros: Faster tree inputs, reduced input space, likely better performance
    - Cons: Can't find the importance of each input (random forest), training the AE could be super expensive and lag out everything
    - This might be a more realistic reach goal than GNNs, lets try this first if there's time. 


### Notes and issues - Tuesday
- It's a simple design problem, but should we have the model be dependent on the dataloader directly, or some direct representation of its data? 
    - Making it dependent on abstract data loaders make it harder to get fine-grained control of what data we're using, so train-test-validate splits are less clear. 
    - The direct meaning of a flat input vector is basically useless to anyone trying to use this package, as the meaning of the vector changes depending on what features you're using for the model.
    - The latter point is- probably -more important in an industry setting, though I don't have much experience here. I'm guessing this is the case: it's more important to have it be readable than to expose fine grain control over training and testing. Therefore: We will add train/test/validate accessors with rng seeds to the dataloader class and implement reasonable defaults. 

- It may be useful to implement predict_pymatgen(pymatgen_structure, extras={}) to facilitate easy prediction with custom structures, but this may be out of scope for now. **Future design will need to consider how exactly novel predictions are going to be made and where from / in what format that data will arrive in.**  

- Clearly, based on parity plots alone we are not predicting well. Though we technically satisfy the "reasonable value" property set forth by the original request (that is, values are in the right ballpark and will interpolate between pure systems), these techniques are a far stretch from being usable in a real environment. We have a number of options and conclusions as of today.
    - Random forest slightly out-performs gradient boosting using these features for large datasets and likely all sizes as well. 
    - Foregoing coulomb matrix methods and attempting more direct- but less informative -features such as simply a list of atomic numbers might assist with these models.
    - Time permitting, an autoencoder would be very interesting to pursue, as it would fix two issues at once: High dimensionality and variable input sizes (necessitating padding).

    - We can check the model for memorization by plotting training parity as well: Thus far memorization does not appear to be occurring with fewer than 100 estimators. This is consistent for smaller datasets. (see below)

    - Since we're padding with zeros and eigenvalues encode structural information, maybe we should ensure the sorting is *descending* and not *ascending*? 
        - Yep! This slightly improves the testing results in both models, almost certainly because larger EVs represent more prevalent informaiton. ]
    
    - I'll go ahead and implement something to find the best hyperparameters and move on from this, as it's a proof of concept. Fixing the descriptor with an autoencoder or another better featurization is likely the core issue, everything else is just slapping band-aids on the problem. 

    - Optimizing hyperparameters definitely doesn't avoid memorization to any degree, and possibly worsens it. 

- Overall conclusion: **Switching from eigenvalues to atomic numbers improves both models RMSE by almost exaclty 0.02, suggesting that the models are struggling to learn to represent these features. As a result, autoencoders present the most viable way to improve these predictions, but implementing them may be out of scope within the time limits of this project**

## Wednesday (~~Finish up model, test model for goals (realistic band gap values, parity, etc)~~ -> Document what exists now, begin understanding more about why the current set of features fail to adequately predict band gap, and begin improving the results.)
### Day plan
- We have a model that overfits/memorizes with larger sizes. Yesterday, I determined that the coulomb matrix eigenvalues, while good descriptors, are too noisy and too ubiquitous to provide the model with generalizable information. Instead, I've come up with a number of options that may avoid this. 
    - Truncating the eigenvalue descriptors to the top N values, possibly cutting out unimportant noise
    - Switching to other (limited) structural information, such as scoring aggregates for the atoms in the lattice + some lattice parameters that are already there
    - Autoencoder representation of the model input (reach, out of time)
    - Graph representation of the model input (very much a reach, but it's here to be a future option)

    ___
    - The first option, truncating the eigenvalues, is the most likely to result in success. The overfitting in the context of the large number of descriptors, as well as the large amount of needed padding, likely happens as the descriptors become "fingerprints" rather than generalizable information. Reducing descriptors therefore is promising.
    

- The classes and function need documentation, namely each model / especially each object that will be exposed in a theoretical API 
- The scripts that provide plots need to be separated from the (rather unorganized) tests, as they're a demonstration rather than tests
- A readme needs to be written to provide instructions on installing and running the PoC



## Thursday (Visualization of results + deal with possible issues or hangups that may have happened Mon-Wed)
### Day plan
- Finish up changed needed to better reduce overfitting
- Document everything in MaterialDataLoader.py
- Implement some more ideas I had
    - Coordination numbers for structural proxies, at various radius
    - Element vectors
    - Dimensionless element vectors


### Notes and issues
- The coulomb matrix is proving to provide no feature importance and is acting more as a fingerprint, as suspected. 
- Instead, we can feed in a list of atomic numbers and some aggregate for connectivity in lieu of a graph
- I think right now we could just get get average coordination numbers using a variety of radii, and also the stdev to tell the model how "spread out" the packing is. 
- Clearly more useful descriptors are needed. The model really cares about mass (not atomic) density. I expected the latter and have not yet found literature to back this. 
- Supplying composition by elemental fractional composition would fix the padding issues 
    - I don't think it would make sense to supply elements that the dataset hasn't previoulsy seen, so we will need to bring attention to the limitation this introduces, lets see if it's worth including. Since chemsys searches are exclusive, they should work better as they have lower atom diversity. 
        - This also needs to be communicated, as if true, it be very important for this models correct usage. 

### Results thus far
- Wow! RMSE 0.44 -> 0.29 with the new descriptors and ways of including structure. This works, but the parity is still indicating overfitting.  
- As expected, limiting the search to a smaller number of atoms (but similarly sized datasets) greatly reduces memorization and RMSE decreases to ~0.08.
- Coulomb eigenvalues significantly harm performance for small datasets as they act like fingerprints, encouraging memorization rather than generalization.
    - CMEs still provide useful information, but more work will need to be done on finding the optimal number of top N eigenvalues to include (hyperparameter tuning problem)
    - I expect the bigger the dataset, the more eigenvalues can be included.

## Friday (last day: wrap up and write comments and instructions.)
### Day plan
- Document PredictorModel.py
- Write an installation readme and point to the demo file. 
- Finish up conclusions in the writeup mdfile. 

