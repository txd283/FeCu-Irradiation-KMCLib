#!/usr/bin/python
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# Used to generate a config.py file for the model
# run script with command './generateConfig.py > config.py'

from KMCLib import *
import numpy

repititons = 5
n = (repititons**3)*2

vacancies = 1
copper = int(round(n*3e-3)) # ~0.3% copper in RPV steel

print "from KMCLib import *"
print("#--------------------------------------------------------------\n")
print("# Number of vacancies (0) = %i"%(vacancies))
print("# Number of Copper (0.1) = %i"%(copper))
print("# Number of Iron (1) = %i\n"%(n - vacancies -copper))
print("#--------------------------------------------------------------\n")
# Define the unit cell.
unit_cell = KMCUnitCell(cell_vectors=numpy.array([[2.87,0.0,0.0],
                                                  [0.0,2.87,0.0],
                                                  [0.0,0.0,2.87]]),
                        basis_points=[[0.0,0.0,0.0],
                                      [0.5,0.5,0.5]])

lattice = KMCLattice(unit_cell=unit_cell,
                     repetitions=(repititons,repititons,repititons),
                     periodic=(True,True,True))

types = ['1']*n

# number of vacancies randomly distributed
for Vacancy in range(vacancies):
    pos = int(numpy.random.rand()*n)
    while (types[pos] == "0"):
        pos = int(numpy.random.rand()*n)
    types[pos] = "0"

# number of vacancies randomly distributed

for Cu in range(copper):
    pos = int(numpy.random.rand()*n)
    while (types[pos] == "0.1"):
        pos = int(numpy.random.rand()*n)
    types[pos] = "0.1"

# Setup the configuration.
configuration = KMCConfiguration(lattice=lattice,
                                 types=types,
                                 possible_types=["0","1","0.1"])

# Use the _script() function to get a script that can generate the configuration.
print configuration._script()
