# here the diffusion processes are described
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *

rate = 1.0

#if diffusion is sucessful, swap the atoms
before_Va_Fe = ['0', '1']
after_Va_Fe = ['1', '0']

before_Va_Cu = ['0', '0.1']
after_Va_Cu = ['0.1', '0']

basis_sites = [0]
second_basis_sites = [1]

# setup available processes. for bcc, there are 8 closest neighbours.

coordinates1 = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
coordinates2 = [[0.0, 0.0, 0.0], [-0.5, 0.5, 0.5]]
coordinates3 = [[0.0, 0.0, 0.0], [0.5, -0.5, 0.5]]
coordinates4 = [[0.0, 0.0, 0.0], [0.5, 0.5, -0.5]]
coordinates5 = [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.5]]
coordinates6 = [[0.0, 0.0, 0.0], [0.5, -0.5, -0.5]]
coordinates7 = [[0.0, 0.0, 0.0], [-0.5, 0.5, -0.5]]
coordinates8 = [[0.0, 0.0, 0.0], [-0.5, -0.5, -0.5]]

# there are 6 second nearest neighbours

coordinates9 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
coordinates10 = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
coordinates11 = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
coordinates12 = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
coordinates13 = [[0.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
coordinates14 = [[0.0, 0.0, 0.0], [0.0, 0.0, -1.0]]

move_vector = None

processes = []

# Va --> Fe diffusion processes for both first and second neighbours

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                             elements_before=before_Va_Fe,
                             elements_after=after_Va_Fe,
                             move_vectors=move_vector,
                             basis_sites=basis_sites,
                             rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

#nearest neighbours, second_basis_sites

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
# add second nearest neighbours in.
# not included as the barrier for second nearest neighbour is very large compared to nearest neighbours.
"""
processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                            
processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
"""
# Va --> Cu diffusion processes for both first

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                             elements_before=before_Va_Cu,
                             elements_after=after_Va_Cu,
                             move_vectors=move_vector,
                             basis_sites=basis_sites,
                             rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=basis_sites,
                            rate_constant=rate))

#nearest neighbours, second_basis_sites

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

"""
processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                            
processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=None,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

"""

# Construct the interactions object.
interactions = KMCInteractions(processes=processes,
                               implicit_wildcards=True)