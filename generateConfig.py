# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import numpy

# Define the unit cell.
unit_cell = KMCUnitCell(cell_vectors=numpy.array([[2.87,0.0,0.0],
                                                  [0.0,2.87,0.0],
                                                  [0.0,0.0,2.87]]),
                        basis_points=[[0.0,0.0,0.0],
                                      [0.5,0.5,0.5]])

lattice = KMCLattice(unit_cell=unit_cell,
                     repetitions=(15,15,15),
                     periodic=(True,True,True))

types = ['Fe']*6750
for i in range(200):
    pos = int(numpy.random.rand()*6750)
    while (types[pos] == "V"):
        pos = int(numpy.random.rand()*6750)
    types[pos] = "V"

# Setup the configuration.
configuration = KMCConfiguration(lattice=lattice,
                                 types=types,
                                 possible_types=["V","Fe"])

# Use the _script() function to get a script that can generate the configuration.
print "from KMCLib import *"
print configuration._script()