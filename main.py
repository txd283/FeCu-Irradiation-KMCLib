# KMCLib application for diffusion of vacancies in pure bcc Fe material.
# Author = Thomas Davis
# email = txd283@bham.ac.uk
# University of Birmingham

from KMCLib import *
import numpy

# add configuration and interactions from files
configuration = KMCConfigurationFromScript("configuration.py")
interactions = KMCInteractionsFromScript("processes.py")

# Generate the KMC model to run.
model = KMCLatticeModel(configuration=configuration,
                        interactions=interactions)

# enable MSD tracking                      
#msd_analysis = OnTheFlyMSD(history_steps=200,
#                           n_bins=100,
#                           t_max=2500.0,
#                           track_type="V")                        
                        
# seed=None uses wall-clock time
control_parameters = KMCControlParameters(number_of_steps=10,
                                          dump_interval=1,
                                          analysis_interval=100,
                                          seed=1234)              

# Run the model and save the atom poisitons to file
model.run(control_parameters=control_parameters,
          trajectory_filename="output.py")
          #trajectory_type = 'lattice',
          #analysis=[msd_analysis])
          
# Save the MSD data to file             
#with open('MSD.data', 'w') as f:
#    msd_analysis.printResults(f)     
          