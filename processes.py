# here the diffusion processes are described
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import math

# need to define a rate, but will be altered by the custom rate class
rate = 1.0

#if diffusion is sucessful, swap the atoms
before_diff = ['V','Fe']
after_diff = ['Fe','V']


basis_sites = [0,1]

# setup available processes. for bcc, there are 8 closest neighbours.

coordinates1=[[0.0,0.0,0.0],[0.5, 0.5, 0.5]]
coordinates2=[[0.0,0.0,0.0],[-0.5, 0.5, 0.5]]
coordinates3=[[0.0,0.0,0.0],[0.5, -0.5, 0.5]]
coordinates4=[[0.0,0.0,0.0],[0.5, 0.5, -0.5]]
coordinates5=[[0.0,0.0,0.0],[-0.5, -0.5, 0.5]]
coordinates6=[[0.0,0.0,0.0],[0.5, -0.5, -0.5]]
coordinates7=[[0.0,0.0,0.0],[-0.5, 0.5, -0.5]]
coordinates8=[[0.0,0.0,0.0],[-0.5, -0.5, -0.5]]

# there are 6 second nearest neighbours

coordinates9=[[0.0,0.0,0.0],[1.0, 0.0, 0.0]]
coordinates10=[[0.0,0.0,0.0],[0.0, 1.0, 0.0]]
coordinates11=[[0.0,0.0,0.0],[0.0, 0.0, 1.0]]
coordinates12=[[0.0,0.0,0.0],[-1.0, 0.0, 0.0]]
coordinates13=[[0.0,0.0,0.0],[0.0, -1.0, 0.0]]
coordinates14=[[0.0,0.0,0.0],[0.0, 0.0, -1.0]]

#move_vectors = [(0,0,0),(1,1,1)]

processes = []

processes.append ( KMCProcess(coordinates=coordinates1,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates2,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates3,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates4,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )

processes.append ( KMCProcess(coordinates=coordinates5,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )

processes.append ( KMCProcess(coordinates=coordinates6,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates7,
                 elements_before=before_diff,
                 elements_after=after_diff,
                 move_vectors=None,
                 basis_sites=basis_sites,
                 rate_constant=rate) )
                 
processes.append ( KMCProcess(coordinates=coordinates8,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
# add second nearest neighbours in

processes.append ( KMCProcess(coordinates=coordinates9,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates10,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )

processes.append ( KMCProcess(coordinates=coordinates11,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )                                                

processes.append ( KMCProcess(coordinates=coordinates12,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )

processes.append ( KMCProcess(coordinates=coordinates13,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
                
processes.append ( KMCProcess(coordinates=coordinates14,
                elements_before=before_diff,
                elements_after=after_diff,
                move_vectors=None,
                basis_sites=basis_sites,
                rate_constant=rate) )
    
# Construct the interactions object.
interactions = KMCInteractions(processes=processes,
                               implicit_wildcards=True)