# here a custom rate is calculated for unique events in the simulation
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
from math import floor
import numpy
import math

# values required for vacancy diffusion, energy in eV, T in K, v is the jump attemps in s^-1

E_m = 0.65
k = 0.862e-4
T = 560
v = 1e13
kT = k*T

# binding energies
E_bVaVa = 0.171
E_bCuCu = 0.235
E_bCuVa = 0.249

# pair interaction values
e_FeVa1 = -0.191
e_FeVa2 = -0.096
e_VaVa1 = 0.255
e_VaVa2 = -0.047
e_CuVa1 = -0.247
e_CuVa2 = -0.206
e_FeCu1 = -0.585
e_FeCu2 = -0.326
e_CuCu1 = -0.627
e_CuCu2 = -0.314


NN1 = [[  1,   2,   3,   4,   5,   6,   7,   8,],
 [  0,   9,  10,  11,  15,  16,  19,  51,],
 [  0,   9,  10,  12,  15,  17,  20,  52,],
 [  0,   9,  11,  13,  16,  18,  21,  53,],
 [  0,   9,  12,  13,  17,  18,  22,  54,],
 [  0,  10,  11,  14,  19,  23,  24,  55,],
 [  0,  10,  12,  14,  20,  23,  25,  56,],
 [  0,  11,  13,  14,  21,  24,  26,  57,],
 [  0,  12,  13,  14,  22,  25,  26,  58,],
 [  1,   2,   3,   4,  27,  28,  29,  30,],
 [  1,   2,   5,   6,  31,  32,  39,  40,],
 [  1,   3,   5,   7,  33,  35,  41,  43,],
 [  2,   4,   6,   8,  34,  36,  42,  44,],
 [  3,   4,   7,   8,  37,  38,  45,  46,],
 [  5,   6,   7,   8,  47,  48,  49,  50,]]
 
 
NN2 = [[  9,  10,  11,  12,  13,  14,],
 [  2,   3,   5,  27,  31,  33,],
 [  1,   4,   6,  28,  32,  34,],
 [  1,   4,   7,  29,  35,  37,],
 [  2,   3,   8,  30,  36,  38,],
 [  1,   6,   7,  39,  41,  47,],
 [  2,   5,   8,  40,  42,  48,],
 [  3,   5,   8,  43,  45,  49,],
 [  4,   6,   7,  44,  46,  50,],
 [  0,  15,  16,  17,  18,  59,],
 [  0,  15,  19,  20,  23,  60,],
 [  0,  16,  19,  21,  24,  61,],
 [  0,  17,  20,  22,  25,  62,],
 [  0,  18,  21,  22,  26,  63,],
 [  0,  23,  24,  25,  26,  64,]]

class CustomRateCalculator(KMCRateCalculatorPlugin):
        
    def initialize(self):
        # used for calculating how many times rate fuction is called.
        self._times_called = 0

    def rate(self, geometry, types_before, types_after, rate_constant, process_number, coordinate):
     
        self._times_called += 1
        print("----------------------- CHECK -----------------------")
        print("Iteration = %i"%(self._times_called))
        
        # find the new position of the moved atom
        for i in range(1,len(types_before)):
            if (float(types_before[i])-float(types_after[i])) != 0.0:
                new_position = i
                break
        
        # types_before first neighbours
        
        b_first_neighbours_Fe = 0.0
        b_first_neighbours_Cu = 0.0
        b_sum_first_FeCu = 0.0
        
        for i in range(8):
            b_sum_first_FeCu += float(types_before[int(NN1[0][i])])
            b_first_neighbours_Fe = floor(b_sum_first_FeCu)
        
        b_first_neighbours_Cu = (b_sum_first_FeCu - floor(b_sum_first_FeCu))*10.0

        b_first_neighbours_Va = abs(8.0 - b_first_neighbours_Fe - b_first_neighbours_Cu)
  
        # types_before second neighbours
  
        b_second_neighbours_Fe = 0.0
        b_second_neighbours_Cu = 0.0
        b_sum_second_FeCu = 0.0
        
        for i in range(6):
            b_sum_second_FeCu += float(types_before[int(NN2[0][i])])
            b_second_neighbours_Fe = floor(b_sum_second_FeCu)
        
        b_second_neighbours_Cu = (b_sum_second_FeCu - floor(b_sum_second_FeCu))*10    
       
        b_second_neighbours_Va = abs(6.0 - b_second_neighbours_Fe - b_second_neighbours_Cu)
        
        # types_after first neighbours
        
        a_first_neighbours_Fe = 0.0
        a_first_neighbours_Cu = 0.0
        a_sum_first_FeCu = 0.0
        
        for i in range(8):
            a_sum_first_FeCu += float(types_after[int(NN1[new_position][i])])
            a_first_neighbours_Fe = floor(a_sum_first_FeCu)
            
        a_first_neighbours_Cu = (a_sum_first_FeCu - floor(a_sum_first_FeCu))*10
        
        a_first_neighbours_Va = abs(8.0 - a_first_neighbours_Fe - a_first_neighbours_Cu)
        
        # types_after second neighbours

        a_second_neighbours_Fe = 0.0
        a_second_neighbours_Cu = 0.0
        a_sum_second_FeCu = 0.0
        
        for i in range(6):
            a_sum_second_FeCu += float(types_after[int(NN2[new_position][i])])
            a_second_neighbours_Fe = floor(a_sum_second_FeCu)
        
        a_second_neighbours_Cu = (a_sum_second_FeCu - floor(a_sum_second_FeCu))*10
         
        a_second_neighbours_Va = abs(6.0 - a_second_neighbours_Fe - a_second_neighbours_Cu)
        
        # energy calculation
       
        E = E_m 
        
        # rate calculation
        
        rate = v*math.exp(-E/kT)
 
        # print out
        
        print("Va moved to %.0f\n"%(new_position))
        print("b_first_neighbours_Va = %.0f"%(b_first_neighbours_Va))
        print("b_second_neighbours_Va = %.0f"%(b_second_neighbours_Va))
        print("b_first_neighbours_Fe = %.0f"%(b_first_neighbours_Fe))
        print("b_second_neighbours_Fe = %.0f"%(b_second_neighbours_Fe))
        print("b_first_neighbours_Cu = %.0f"%(b_first_neighbours_Cu))
        print("b_second_neighbours_Cu = %.0f\n"%(b_second_neighbours_Cu))
        
        print("a_first_neighbours_Va = %.0f"%(a_first_neighbours_Va))
        print("a_second_neighbours_Va = %.0f"%(a_second_neighbours_Va))
        print("a_first_neighbours_Fe = %.0f"%(a_first_neighbours_Fe))
        print("a_second_neighbours_Fe = %.0f"%(a_second_neighbours_Fe))
        print("a_first_neighbours_Cu = %.0f"%(a_first_neighbours_Cu))
        print("a_second_neighbours_Cu = %.0f\n"%(a_second_neighbours_Cu))
        
        print("E = %.4f eV"%(E))        
        print ("Rate = %.3f\n"%(rate))   
        
        # return the new rate value
        return rate
        
    def cutoff(self):
        # cutoff value for types_before and types_after lattice points. 2.0 = two supercells.
        return 2.0