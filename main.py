# KMCLib application for diffusion of vacancies in pure bcc Fe material.
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# to run the simulation, go into command line and type 'python main.py'

#from mpi4py import MPI
from KMCLib import *
from customRateCalculator import *

# add configuration and interactions from files
configuration = KMCConfigurationFromScript("config.py")
interactions = KMCInteractionsFromScript("processes.py")


# Set the rate calculator which includes vacancy clustering
interactions.setRateCalculator(rate_calculator=CustomRateCalculator)


# Generate the KMC model to run.
model = KMCLatticeModel(configuration=configuration,
                        interactions=interactions)


# set number_of_steps to wanted value
# dump_interval set to a value to dump the current locations.                                          
# seed=None uses wall-clock time
control_parameters = KMCControlParameters(number_of_steps=10000,
                                          dump_interval=10000,
                                          seed=None)              

# here define type of trajectory you want to dump. I have created a script which converts the results1.py to CFG files for visulatiation in atomeye by convertToCFG.py in folder UsefulCode
# Run the model and save the atom poisitons to file
model.run(control_parameters=control_parameters,
              trajectory_filename="results/results1.py",
              trajectory_type = 'lattice')
              
