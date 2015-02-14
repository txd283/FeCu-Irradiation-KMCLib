# here the configuration for the lattice will be determined 
# Author = Thomas Davis
# email = txd283@bham.ac.uk
# University of Birmingham

from KMCLib import *
import numpy
import math
# Define a lattice - Fe lattice with a = 2.87 A
# basis_points there is an atom in the middle of the unit cell at [2,2,2]

unit_cell = KMCUnitCell(cell_vectors=numpy.array([[2.87,0.0,0.0],
                                                  [0.0,2.87,0.0],
                                                  [0.0,0.0,2.87]]),
                        basis_points=[[0.0,0.0,0.0],
                                      [0.5,0.5,0.5]])

lattice = KMCLattice(unit_cell=unit_cell,
                     repetitions=(10,10,10),
                     periodic=(True,True,True))


# Types of atoms. There are total 2 atoms in a unit cell of bcc => 10x10x10 = 1,000 unit cells = 2000 atoms
types = ["Fe"]*2000

# Swap parts of the array with vacancies, V.
types[20] = "V"
types[22] = "V"
types[19] = "V"
types[111] = "V"
types[44] = "V"

# Generate the configuration. With the lattice, types and possible types.
config = KMCConfiguration(lattice=lattice,
                          types=types,
                          possible_types=["Fe","V"])