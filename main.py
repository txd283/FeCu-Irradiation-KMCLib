# KMCLib application for diffusion of vacancies in pure bcc Fe material.
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
from customRates import *
import numpy
import os

# add configuration and interactions from files
configuration = KMCConfigurationFromScript("config.py")
interactions = KMCInteractionsFromScript("processes.py")

# Set the rate calculator which includes vacancy clustering
interactions.setRateCalculator(rate_calculator=CustomRateCalculator)

# Generate the KMC model to run.
model = KMCLatticeModel(configuration=configuration,
                        interactions=interactions)

# enable MSD tracking                      
msd_analysis = OnTheFlyMSD(history_steps=1000,
                           n_bins=100,
                           t_max=2500.0,
                           track_type="V")                        
                        
# seed=None uses wall-clock time
control_parameters = KMCControlParameters(number_of_steps=10,
                                          dump_interval=1,
                                          analysis_interval=1,
                                          seed=111)              

# Run the model and save the atom poisitons to file
model.run(control_parameters=control_parameters,
          trajectory_filename="results/all.cfg",
          trajectory_type = 'cfg',
          analysis=[msd_analysis])
          
# Save the MSD data to file             
with open('results/MSD.data', 'w') as f:
    msd_analysis.printResults(f)     

#split the cfg into different timesteps for atomeye
numberOfAtoms = 250
# 13 is the line number of the header printed.
lines_per_file = 13 + numberOfAtoms 
number = 0
smallfile = 0
with open('results/all.cfg') as bigfile:
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if smallfile:
                smallfile.close()
            small_filename = 'results/00{}.cfg'.format(number)
            smallfile = open(small_filename, "w")
            number = number + 1
        smallfile.write(line)
    if smallfile:
        smallfile.close()

#open up atomeye with the cfg files generated
#os.system(" open /usr/local/bin/atomeye results/cfg000.cfg")
