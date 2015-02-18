#!/usr/bin/python
# KMCLib application for diffusion of vacancies in pure bcc Fe material.
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
from customRates import *
import numpy

# add configuration and interactions from files
configuration = KMCConfigurationFromScript("config.py")
interactions = KMCInteractionsFromScript("processes.py")


# Set the rate calculator which includes vacancy clustering
interactions.setRateCalculator(rate_calculator=CustomRateCalculator)


# Generate the KMC model to run.
model = KMCLatticeModel(configuration=configuration,
                        interactions=interactions)
                                          
# seed=None uses wall-clock time
control_parameters = KMCControlParameters(number_of_steps=1000000,
                                          dump_interval=20000,
                                          seed=None)              

# Run the model and save the atom poisitons to file
model.run(control_parameters=control_parameters,
          trajectory_filename="results/lattice3.py",
          trajectory_type = 'lattice')
