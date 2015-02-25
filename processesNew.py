# here the diffusion processes are described
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import math

rate = 1.0

#if diffusion is sucessful, swap the atoms
before_Va_Fe = ['0', '1']
after_Va_Fe = ['1', '0']

#before_Va_Cu = ['Va', 'Cu']
#after_Va_Cu = ['Cu', 'Va']

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

move_vector1 = [(0, [0.5, 0.5, 0.5]),
                (1, [-0.5, -0.5, -0.5])]
                
move_vector2 = [(0, [-0.5, 0.5, 0.5]),
                (1, [0.5, -0.5, -0.5])]
                
move_vector3 = [(0, [0.5, -0.5, 0.5]),
                (1, [-0.5, 0.5, -0.5])]

move_vector4 = [(0, [0.5, 0.5, -0.5]),
                (1, [-0.5, -0.5, 0.5])]
                
move_vector5 = [(0, [-0.5, -0.5, 0.5]),
                (1, [0.5, 0.5, -0.5])]
                
move_vector6 = [(0, [0.5, -0.5, -0.5]),
                (1, [-0.5, 0.5, 0.5])]
                
move_vector7 = [(0, [-0.5, 0.5, -0.5]),
                (1, [0.5, -0.5, 0.5])]
                
move_vector8 = [(0, [-0.5, -0.5, -0.5]),
                (1, [0.5, 0.5, 0.5])]

move_vector9 = [(0, [1.0, 0.0, 0.0]),
                (1, [-1.0, 0.0, 0.0])]
                
move_vector10 = [(0, [0.0, 1.0, 0.0]),
                 (1, [0.0, -1.0, 0.0])]
                 
move_vector11 = [(0, [0.0, 0.0, 1.0]),
                 (1, [0.0, 0.0, -1.0])]
                 
move_vector12 = [(0, [-1.0, 0.0, 0.0]),
                 (1, [1.0, 0.0, 0.0])]
                 
move_vector13 = [(0, [0.0, -1.0, 0.0]),
                 (1, [0.0, 1.0, 0.0])]
                 
move_vector14 = [(0, [0.0, 0.0, -1.0]),
                 (1, [0.0, 0.0, -1.0])]

move_vector_second_1 = [(1, (0.5, 0.5, 0.5))]
move_vector_second_2 = [(1, (-0.5, 0.5, 0.5))]
move_vector_second_3 = [(1, (0.5, -0.5, 0.5))]
move_vector_second_4 = [(1, (0.5, 0.5, -0.5))]
move_vector_second_5 = [(1, (-0.5, -0.5, 0.5))]
move_vector_second_6 = [(1, (0.5, -0.5, -0.5))]
move_vector_second_7 = [(1, (-0.5, 0.5, -0.5))]
move_vector_second_8 = [(1, (-0.5, -0.5, -0.5))]

move_vector_second_9 = [(1, (1.0, 0.0, 0.0))]
move_vector_second_10 = [(1, (0.0, 1.0, 0.0))]
move_vector_second_11 = [(1, (0.0, 0.0, 1.0))]
move_vector_second_12 = [(1, (-1.0, 0.0, 0.0))]
move_vector_second_13 = [(1, (0.0, -1.0, 0.0))]
move_vector_second_14 = [(1, (0.0, 0.0, -1.0))]

processes = []

# Va --> Fe diffusion processes for both first and second neighbours

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector1,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector2,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector3,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector4,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector5,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector6,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                             elements_before=before_Va_Fe,
                             elements_after=after_Va_Fe,
                             move_vectors=move_vector7,
                             basis_sites=basis_sites,
                             rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector8,
                            basis_sites=basis_sites,
                            rate_constant=rate))

#nearest neighbours, second_basis_sites
""""
processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_1,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_2,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_3,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_4,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_5,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_6,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_7,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_8,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
"""           
# add second nearest neighbours in

processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector9,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector10,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector11,
                            basis_sites=basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector12,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector13,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector14,
                            basis_sites=basis_sites,
                            rate_constant=rate))
""""                            
processes.append(KMCProcess(coordinates=coordinates9,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_9,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates10,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_10,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates11,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_11,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))  

processes.append(KMCProcess(coordinates=coordinates12,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_12,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates13,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_13,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates14,
                            elements_before=before_Va_Fe,
                            elements_after=after_Va_Fe,
                            move_vectors=move_vector_second_14,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
"""
# Va --> Cu diffusion processes for both first and second neighbours
"""
processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector1,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector2,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector3,
                            basis_sites=basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector4,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector5,
                            basis_sites=basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector6,
                            basis_sites=basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                             elements_before=before_Va_Cu,
                             elements_after=after_Va_Cu,
                             move_vectors=move_vector7,
                             basis_sites=basis_sites,
                             rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector8,
                            basis_sites=basis_sites,
                            rate_constant=rate))

#nearest neighbours, second_basis_sites

processes.append(KMCProcess(coordinates=coordinates1,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector1,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates2,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector2,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates3,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector3,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
                
processes.append(KMCProcess(coordinates=coordinates4,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector4,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates5,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector6,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))

processes.append(KMCProcess(coordinates=coordinates6,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector6,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
              
processes.append(KMCProcess(coordinates=coordinates7,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector7,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))
             
processes.append(KMCProcess(coordinates=coordinates8,
                            elements_before=before_Va_Cu,
                            elements_after=after_Va_Cu,
                            move_vectors=move_vector8,
                            basis_sites=second_basis_sites,
                            rate_constant=rate))                        
"""
# Construct the interactions object.
interactions = KMCInteractions(processes=processes,
                               implicit_wildcards=True)