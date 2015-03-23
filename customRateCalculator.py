# Custom rate is calculated for unique events in the simulation
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# This rate calculator has a full implementation of the pair interaction method
# for up to second nearest neighbours
# Cu and Vacancies clustering works
# Only first neighbour hops are included in the processes.py, and thus no code 
# here has been incoperated for second neighbour hops

from KMCLib import *
from math import floor
import numpy
import math

# values required for vacancy diffusion, energy in eV, T in K, v is the jump in s^-1

E_m_Fe = 0.722
E_m_Cu = 0.50
k = 0.862e-4
T = 563
v_Fe = 9.79e12
v_Cu = 7.16e12
kT = k*T

# pair interaction values
e_FeFe1 = -0.778
e_FeFe2 = -0.389
e_VaFe1 = -0.191
e_VaFe2 = -0.096
e_VaVa1 = 0.225
e_VaVa2 = -0.047
e_VaCu1 = -0.247
e_VaCu2 = -0.206
e_FeCu1 = -0.585
e_FeCu2 = -0.326
e_CuCu1 = -0.627
e_CuCu2 = -0.314

# The first nearest neighbours for all atoms in the lattice in types_before and types_after. 
# Used to find local configurations for energy calculations
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
 
# The second nearest neighbours for all atoms in the lattice in types_before and types_after. 
# Used to find local configurations for energy calculations
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
    
   # uncomment if you want a counter for the number of times the fuction rate() is called -- for diagnositics.            
   # def initialize(self):
        # used for calculating how many times rate fuction is called.
        #self._times_called = 0

    def rate(self, geometry, types_before, types_after, rate_constant, process_number, coordinate):
        # see above -- diagnositics
        #self._times_called += 1
                
        # find the new position of the moved atom
        for i in range(1,len(types_before)):
            if (float(types_before[i])-float(types_after[i])) != 0.0:
                new_position = i
                break
        
        # define variables for the pair interaction. N_FeFe1_b is the number of Fe-Fe bonds before the move and 'a' stand for after the move.
        
        N_FeFe1_b = 0.0
        N_FeFe1_a = 0.0
        N_FeFe2_b = 0.0
        N_FeFe2_a = 0.0
        
        N_CuCu1_b = 0.0
        N_CuCu1_a = 0.0
        N_CuCu2_b = 0.0
        N_CuCu2_a = 0.0
        
        N_VaVa1_b = 0.0
        N_VaVa1_a = 0.0
        N_VaVa2_b = 0.0
        N_VaVa2_a = 0.0
        
        N_VaFe1_b = 0.0
        N_VaFe1_a = 0.0
        N_VaFe2_b = 0.0
        N_VaFe2_a = 0.0
        
        N_FeVa1_b = 0.0
        N_FeVa1_a = 0.0
        N_FeVa2_b = 0.0
        N_FeVa2_a = 0.0
        
        N_VaCu1_b = 0.0
        N_VaCu1_a = 0.0
        N_VaCu2_b = 0.0
        N_VaCu2_a = 0.0
          
        N_CuVa1_b = 0.0
        N_CuVa1_a = 0.0
        N_CuVa2_b = 0.0
        N_CuVa2_a = 0.0      
        
        N_FeCu1_b = 0.0
        N_FeCu1_a = 0.0
        N_FeCu2_b = 0.0
        N_FeCu2_a = 0.0
        
        N_CuFe1_b = 0.0
        N_CuFe1_a = 0.0
        N_CuFe2_b = 0.0
        N_CuFe2_a = 0.0
        
        # find first neighbours of Va before move
        count = 0.0
        # count the number of bonds at position 0 at all the possible nearest neighbours. Uses NN1 array.
        for i in range(8):
            count += float(types_before[int(NN1[0][i])])
        # floor will reveal the number of  1st nearest neighbour Va-Fe bonds
        N_VaFe1_b = floor(count)
        # will reveal the number of 1st nearest neighbour Va-Cu bonds
        N_VaCu1_b = (count - floor(count))*10.0
        # remaining values will be  1st nearest neighbour Va-Va bonds
        N_VaVa1_b = abs(8.0 - N_VaFe1_b - N_VaCu1_b)

        # same method above, but now for 2nd nearest neighbours. uses NN2 array.
        count = 0.0
        for i in range(6):
            count += float(types_before[int(NN2[0][i])])
        N_VaFe2_b = floor(count)
        N_VaCu2_b = (count - floor(count))*10.0
        N_VaVa2_b = abs(6.0 - N_VaFe2_b - N_VaCu2_b)
        
        # find first neighbours of Va after move
        count = 0.0
        for i in range(8):
            count += float(types_after[int(NN1[new_position][i])])
        N_VaFe1_a = floor(count)
        N_VaCu1_a = (count - floor(count))*10.0
        N_VaVa1_a = abs(8.0 - N_VaFe1_a - N_VaCu1_a)
  
        # find second neighbours of Va after move
        count = 0.0
        for i in range(6):
            count += float(types_after[int(NN2[new_position][i])])
        N_VaFe2_a = floor(count)
        N_VaCu2_a = (count - floor(count))*10.0
        N_VaVa2_a = abs(6.0 - N_VaFe2_a - N_VaCu2_a)
        
        # Find what atom the Va is swapping with - either a Fe (1) or Cu(0.1)
        if types_after[0] == "1":
            
            # find first neighbours of Fe before move
            count = 0.0
            for i in range(8):
                count += float(types_before[int(NN1[new_position][i])])
            N_FeFe1_b = floor(count)
            N_FeCu1_b = (count - floor(count))*10.0
            N_FeVa1_b = abs(8.0 - N_FeFe1_b - N_FeCu1_b)
            
            # find second neighbours of Fe before move
            count = 0.0
            for i in range(6):
                count += float(types_before[int(NN2[new_position][i])])
            N_FeFe2_b = floor(count)
            N_FeCu2_b = (count - floor(count))*10.0
            N_FeVa2_b = abs(6.0 - N_FeFe2_b - N_FeCu2_b)
        
            # find first neighbours of Fe after move
            count = 0.0
            for i in range(8):
                count += float(types_after[int(NN1[0][i])])
            N_FeFe1_a = floor(count)
            N_FeCu1_a = (count - floor(count))*10.0
            N_FeVa1_a = abs(8.0 - N_FeFe1_a - N_FeCu1_a)
  
            # find second neighbours of Fe after move
            count = 0.0
            for i in range(6):
                count += float(types_after[int(NN2[0][i])])
            N_FeFe2_a = floor(count)
            N_FeCu2_a = (count - floor(count))*10.0
            N_FeVa2_a = abs(6.0 - N_FeFe2_a - N_FeCu2_a)
        
        else:
            
            # find first neighbours of Cu before move
            count = 0.0
            for i in range(8):
                count += float(types_before[int(NN1[new_position][i])])
            N_CuFe1_b = floor(count)
            N_CuCu1_b = (count - floor(count))*10.0
            N_CuVa1_b = abs(8.0 - N_CuFe1_b - N_CuCu1_b)
            
            # find second neighbours of Cu before move
            count = 0.0
            for i in range(6):
                count += float(types_before[int(NN2[new_position][i])])
            N_CuFe2_b = floor(count)
            N_CuCu2_b = (count - floor(count))*10.0
            N_CuVa2_b = abs(6.0 - N_CuFe2_b - N_CuCu2_b)
        
            # find first neighbours of Cu after move
            count = 0.0
            for i in range(8):
                count += float(types_after[int(NN1[0][i])])
            N_CuFe1_a = floor(count)
            N_CuCu1_a = (count - floor(count))*10.0
            N_CuVa1_a = abs(8.0 - N_CuFe1_a - N_CuCu1_a)
  
            # find second neighbours of Cu after move
            count = 0.0
            for i in range(6):
                count += float(types_after[int(NN2[0][i])])
            N_CuFe2_a = floor(count)
            N_CuCu2_a = (count - floor(count))*10.0
            N_CuVa2_a = abs(6.0 - N_CuFe2_a - N_CuCu2_a)
        
        
        # find the difference before and after the jump bonds.
        D_N_FeFe1 = N_FeFe1_a - N_FeFe1_b
        D_N_FeFe2 = N_FeFe2_a - N_FeFe2_b
        D_N_CuCu1 = N_CuCu1_a - N_CuCu1_b
        D_N_CuCu2 = N_CuCu2_a - N_CuCu2_b
        D_N_VaVa1 = N_VaVa1_a - N_VaVa1_b
        D_N_VaVa2 = N_VaVa2_a - N_VaVa2_b
        D_N_VaFe1 = N_VaFe1_a + N_FeVa1_a - N_VaFe1_b - N_FeVa1_b
        D_N_VaFe2 = N_VaFe2_a + N_FeVa2_a - N_VaFe2_b - N_FeVa2_b    
        D_N_VaCu1 = N_VaCu1_a + N_CuVa1_a - N_VaCu1_b - N_CuVa1_b
        D_N_VaCu2 = N_VaCu2_a + N_CuVa2_a - N_VaCu2_b - N_CuVa2_b
        D_N_FeCu1 = N_FeCu1_a + N_CuFe1_a - N_FeCu1_b - N_CuFe1_b
        D_N_FeCu2 = N_FeCu2_a + N_CuFe2_a - N_FeCu2_b - N_CuFe2_b
        
        # binding energy calculation
        E_b = (D_N_VaVa1*e_VaVa1 +
               D_N_VaFe1*e_VaFe1 + 
               D_N_FeFe1*e_FeFe1 + 
               D_N_VaVa2*e_VaVa2 + 
               D_N_VaFe2*e_VaFe2 + 
               D_N_FeFe2*e_FeFe2 + 
               D_N_VaCu1*e_VaCu1 + 
               D_N_VaCu2*e_VaCu2 + 
               D_N_FeCu1*e_FeCu1 + 
               D_N_FeCu2*e_FeCu2 + 
               D_N_CuCu1*e_CuCu1 + 
               D_N_CuCu2*e_CuCu2)/2
    
        # if the atom is eith Fe (1) or Cu (0.1) and calculate the rate
        if types_after[0] == '1':
            E = E_m_Fe + E_b
            rate = v_Fe*math.exp(-E/kT)
            
        else:
            E = E_m_Cu + E_b
            rate = v_Cu*math.exp(-E/kT)
        
        # commented out this code -- it is the implementation of second nearest neighbours hops. need to define new E_m and v_Fe/Cu values, as they are different for first and second neighbour hops.  
        """
        distance_sq = (geometry[new_position][0] - geometry[0][0] )**2 + (geometry[new_position][1] - geometry[0][1] )**2 + (geometry[new_position][2] - geometry[0][2] )**2        
    
        if types_after[0] == '1' and distance_sq == 0.75: 
            E = E_m_Fe_1 + E_b
            rate = v_Fe_1*math.exp(-E/kT)
            print("types_after[0] == '1' and distance_sq == 0.75:")
        elif types_after[0] == '1' and distance_sq == 1.0:
            E = E_m_Fe_2 + E_b
            rate = v_Fe_2*math.exp(-E/kT)
            print("types_after[0] == '1' and distance_sq == 1.0:")
            
        elif types_after[0] == '0.1' and distance_sq == 0.75:
            E = E_m_Cu_1 + E_b
            rate = v_Cu_1*math.exp(-E/kT)
            print("types_after[0] == '0.1' and distance_sq == 0.75:")
        elif types_after[0] == '0.1' and distance_sq == 1.0:
            E = E_m_Cu_2 + E_b
            rate = v_Cu_2*math.exp(-E/kT)
            print("types_after[0] == '0.1' and distance_sq == 1.0:")
            
        
        # print out useful variables for diagnostics 
        print("----------------------- CHECK -----------------------")
        print("Iteration = %i"%(self._times_called))
        print("Vacancy moved to %i\n"%(new_position))
        print("D_N_FeFe1 = %.0f"%(D_N_FeFe1))
        print("D_N_VaFe1 = %.0f"%(D_N_VaFe1))
        print("D_N_VaVa1 = %.0f"%(D_N_VaVa1))
        print("D_N_VaCu1 = %.0f"%(D_N_VaCu1))
        print("D_N_FeCu1 = %.0f"%(D_N_FeCu1))
        print("D_N_CuCu1 = %.0f\n"%(D_N_CuCu1))
        print("D_N_FeFe2 = %.0f"%(D_N_FeFe2))
        print("D_N_VaFe2 = %.0f"%(D_N_VaFe2))
        print("D_N_VaVa2 = %.0f"%(D_N_VaVa2))
        print("D_N_VaCu2 = %.0f"%(D_N_VaCu2))
        print("D_N_FeCu2 = %.0f"%(D_N_FeCu2))
        print("D_N_CuCu2 = %.0f\n"%(D_N_CuCu2))
        print("E_b = %.2f eV"%(E_b))        
        print("E = %.4f eV"%(E))        
        print ("Rate = %f\n"%(rate))   
        """
        
        # return the new rate value
        return rate
        
    def cutoff(self):
        # cutoff value for types_before and types_after lattice points. 2.0 = two supercells.
        return 2.0