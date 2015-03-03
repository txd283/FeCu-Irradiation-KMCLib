#!/usr/bin/python
# KMCLib application for diffusion of vacancies in pure bcc Fe material.
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
from customRates_VaCu import *
from timer import Timer

# add configuration and interactions from files
configuration = KMCConfigurationFromScript("config.py")
interactions = KMCInteractionsFromScript("processes.py")


# Set the rate calculator which includes vacancy clustering
interactions.setRateCalculator(rate_calculator=CustomRateCalculator)


# Generate the KMC model to run.
model = KMCLatticeModel(configuration=configuration,
                        interactions=interactions)
                                          
# seed=None uses wall-clock time
control_parameters = KMCControlParameters(number_of_steps=0,
                                          dump_interval=1,
                                          seed=None)              

# Run the model and save the atom poisitons to file
with Timer() as t:
    model.run(control_parameters=control_parameters,
              trajectory_filename="results/test.py",
              trajectory_type = 'lattice')
print("Total CPU Time = %s s"%(t.secs))

#print(CustomRateCalculator.initialize)
