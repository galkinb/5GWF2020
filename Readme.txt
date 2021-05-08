This is the source code for the paper "Intelligent Base Station Association for UAV Cellular Users: A Supervised Learning Approach" 
by B. Galkin, R. Amer, E. Fonseca, and L.A. DaSilva, published in the IEEE 5G World Forum, September 2020. The code consists of several files

1. NNTrainingDataset.R generates the dataset that will be used to train the neural network
2. NNTraining.R uses this dataset to train the neural network.
3. NNEvaluation.R compares the performance of our trained neural network against benchmarks across different heights
(the comparison for different densities and beamwidths is a very simple modification to the code, and is ommitted here)

ClosedFormFunctions.R includes the stochastic geometry expressions that are used to generate the nearest BS association benchmark results.
MiscFunctions.R includes miscellaneous functions that are used in the code.