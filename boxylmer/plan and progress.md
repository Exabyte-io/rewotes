# Minimum Goals
- Create a model which will automatically fit and predict material band-gaps given an input. 

# The Plan
## Monday (primarily planning, familiarize with literature and MPRester)
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
- Identify appropriate structure for models and fill out the BandGapPredictor
    - in this limited case the input will be flat, but a 3D structure such as graphs will be more complicated and likely require a different representation or accessor method
    - The loader should hold onto the structure objects, so graph construction is still possible without weird methods / a class rewrite.
- Models: RandomForest (doable), GradientBoosting (doable), Graph Convolution (reach, likely won't have enough time this week)
- Padding -> Decide how to allocate it, how will the max size be determined and stored?
- Possibly, we could use an autoencoder to get over possible padding issues, if the tree based methods don't like the zero padding. 
    - Pros: Faster tree inputs, reduced input space, likely better performance
    - Cons: Can't find the importance of each input (random forest), training the AE could be super expensive and lag out everything
    - This might be a more realistic reach goal than GNNs, lets try this first if there's time.   

## Wednesday (Finish up model, test model for goals (realistic band gap values, parity, etc))

## Thursday (Visualization of results + deal with possible issues or hangups that may have happened Mon-Wed)

## Friday (last day: wrap up and write comments and instructions.)