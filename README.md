# SignalingMatingTypes
Source code for Hadjivasiliou and Pomiankowski, "Evolution of asymmetric gamete signaling and suppressed recombination at the mating type locus". eLife, 2019. http://dx.doi.org/10.7554/eLife.48239

# C++ code
Simulations can be run by compiling all files ending in cpp. main.cpp contains the core of the simulations intiating matrices and a loop with mutation and mating. Signal.cpp contains all functions used to implement signaling between and within cells, mutation and mating. Functions.cpp contains a collection of general functions that are used throughout the code. InitiationFunctions.cpp includes functions used to initialize matrices and WriteFiles.cpp includes functions that produce the output at the end of the simulation.


# Mathematica notebooks
Mathematica was used to analyse simulated data. 

CoEvo_s_r reads in and plots the evolution trajectories for s and r.

Evo_Averaged reads in the evolution trajectories for specified parameter sets and plots the average steady state value for s or r
