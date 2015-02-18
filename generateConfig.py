# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import numpy
import sys

repititons = 40

# Define the unit cell.
unit_cell = KMCUnitCell(cell_vectors=numpy.array([[2.87,0.0,0.0],
                                                  [0.0,2.87,0.0],
                                                  [0.0,0.0,2.87]]),
                        basis_points=[[0.0,0.0,0.0],
                                      [0.5,0.5,0.5]])

lattice = KMCLattice(unit_cell=unit_cell,
                     repetitions=(repititons,repititons,repititons),
                     periodic=(True,True,True))
n = (repititons**3)*2

#print("Generating lattice with %i of atoms..."%(n))
types = ['Fe']*n
for i in range(500):
    pos = int(numpy.random.rand()*n)
    while (types[pos] == "V"):
        pos = int(numpy.random.rand()*n)
    types[pos] = "V"
#print("Done")

# Setup the configuration.
configuration = KMCConfiguration(lattice=lattice,
                                 types=types,
                                 possible_types=["V","Fe"])
# Use the _script() function to get a script that can generate the configuration.
print "from KMCLib import *"
print configuration._script()
#print configuration._atkScript(types)