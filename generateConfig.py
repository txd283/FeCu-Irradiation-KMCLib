#!/usr/bin/python
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# Used to generate a config.py file for the model
# run script with command './generateConfig.py > config.py'

from KMCLib import *
import numpy

# here enter the repititons of the cell.
repititons = 10
n = (repititons**3)*2

# here enter the number of vacancies and copper atoms you would like -- they will be dispersed randomly. Copper could be added with a set number based on the wt%, or override this with an integer value.

vacancies = 1
copper = int(round(n*5e-4)) # 0.05% copper in RPV steel

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

# Here Fe = 1, Vacancies = 0 and Cu = 0.1 for making the custom rate calculator efficient.
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
