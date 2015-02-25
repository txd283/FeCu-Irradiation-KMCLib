#!/usr/bin/python
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import numpy
import sys

repititons = 10
n = (repititons**3)*2

vacancies = 100
copper = int(round(n*3e-3))

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
"""
for Cu in range(copper):
    pos = int(numpy.random.rand()*n)
    while (types[pos] == "Cu"):
        pos = int(numpy.random.rand()*n)
    types[pos] = "Cu"
"""
# Setup the configuration.
configuration = KMCConfiguration(lattice=lattice,
                                 types=types,
                                 possible_types=["0","1","0.1"])

# Use the _script() function to get a script that can generate the configuration.
print "from KMCLib import *"
print configuration._script()
