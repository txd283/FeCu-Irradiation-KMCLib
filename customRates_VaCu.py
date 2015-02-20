# here a custom rate is calculated for unique events in the simulation
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
from interactionCount import counti
import numpy
import math

#### values required for vacancy diffusion, energy in eV, T in K, v is the jump attemps in s^-1
E_m = 0.55
E_b = 0.2
k = 0.862e-4
T = 560
v = 1e13
kT = k*T
diffusion_rate_FeVa = v*math.exp(-E_m/(kT))

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

class CustomRateCalculator(KMCRateCalculatorPlugin):

    def rate(self, geometry, types_before, types_after, rate_constant, process_number, coordinate):

        # 
        # BEFORE EVENT
        #
        # For every movable atom, check the neighbours around it to see if they are V.
        b_first_neighbours_Va = len([ i for i in [types_before[1], types_before[2], types_before[3], types_before[4], types_before[5],types_before[6], types_before[7], types_before[8]] if i == "Va"])
        b_second_neighbours_Va = len([ i for i in [types_before[9], types_before[10], types_before[11], types_before[12], types_before[13],types_before[14]] if i == "Va"])
        
        # For every movable atom, check the neighbours around it to see if they are Cu.
        b_first_neighbours_Cu = len([ i for i in [types_before[1], types_before[2], types_before[3], types_before[4], types_before[5],types_before[6], types_before[7], types_before[8]] if i == "Cu"])
        b_second_neighbours_Cu = len([ i for i in [types_before[9], types_before[10], types_before[11], types_before[12], types_before[13],types_before[14]] if i == "Cu"])
        
        # For every movable atom, check the neighbours around it to see if they are Fe.
        b_first_neighbours_Fe = len([ i for i in [types_before[1], types_before[2], types_before[3], types_before[4], types_before[5],types_before[6], types_before[7], types_before[8]] if i == "Fe"])
        b_second_neighbours_Fe = len([ i for i in [types_before[9], types_before[10], types_before[11], types_before[12], types_before[13],types_before[14]] if i == "Fe"])
        
        #
        # AFTER EVENT
        #
        # For every movable atom, check the neighbours around it to see if they are Va.
        a_first_neighbours_Va = len([ i for i in [types_after[1], types_after[2], types_after[3], types_after[4], types_after[5],types_after[6], types_after[7], types_after[8]] if i == "Va"])
        a_second_neighbours_Va = len([ i for i in [types_after[9], types_after[10], types_after[11], types_after[12], types_after[13],types_after[14]] if i == "Va"])
        
        # For every movable atom, check the neighbours around it to see if they are Cu.
        a_first_neighbours_Cu = len([ i for i in [types_after[1], types_after[2], types_after[3], types_after[4], types_after[5],types_after[6], types_after[7], types_after[8]] if i == "Cu"])
        a_second_neighbours_Cu = len([ i for i in [types_after[9], types_after[10], types_after[11], types_after[12], types_after[13],types_after[14]] if i == "Cu"])
        
        # For every movable atom, check the neighbours around it to see if they are Fe.
        a_first_neighbours_Fe = len([ i for i in [types_after[1], types_after[2], types_after[3], types_after[4], types_after[5],types_after[6], types_after[7], types_after[8]] if i == "Fe"])
        a_second_neighbours_Fe = len([ i for i in [types_after[9], types_after[10], types_after[11], types_after[12], types_after[13],types_after[14]] if i == "Fe"])
        
        print("---------- Check -----------")
        print("b_first_neighbours_Va = %i"%(b_first_neighbours_Va))
        print("b_second_neighbours_Va = %i"%(b_second_neighbours_Va))
        print("b_first_neighbours_Cu = %i"%(b_first_neighbours_Cu))
        print("b_second_neighbours_Cu = %i"%(b_second_neighbours_Cu))
        print("b_first_neighbours_Fe = %i"%(b_first_neighbours_Fe))
        print("b_second_neighbours_Fe = %i\n"%(b_second_neighbours_Fe))
        
        print("a_first_neighbours_Va = %i"%(a_first_neighbours_Va))
        print("a_second_neighbours_Va = %i"%(a_second_neighbours_Va))
        print("a_first_neighbours_Cu = %i"%(a_first_neighbours_Cu))
        print("a_second_neighbours_Cu = %i"%(a_second_neighbours_Cu))
        print("a_first_neighbours_Fe = %i"%(a_first_neighbours_Fe))
        print("a_second_neighbours_Fe = %i\n"%(a_second_neighbours_Fe))
        
        # pair interaction calculation
        
        E_before = b_first_neighbours_Fe*e_FeVa1 + b_second_neighbours_Fe*e_FeVa2 + b_first_neighbours_Va*e_VaVa1 + b_second_neighbours_Va*e_VaVa2 + b_second_neighbours_Cu*e_CuVa1 + b_second_neighbours_Cu*e_CuVa2 + b_second_neighbours_Cu*e_FeCu1 + b_second_neighbours_Cu*e_FeCu2
        
        E_after = a_first_neighbours_Fe*e_FeVa1 + a_second_neighbours_Fe*e_FeVa2 + a_first_neighbours_Va*e_VaVa1 + a_second_neighbours_Va*e_VaVa2 + a_second_neighbours_Cu*e_CuVa1 + a_second_neighbours_Cu*e_CuVa2 + a_second_neighbours_Cu*e_FeCu1 + a_second_neighbours_Cu*e_FeCu2
        
        E_pi = 0.5*(E_after - E_before)
        
        print("E_before = %f"%(E_before))
        print("E_after = %f"%(E_after))
        print("E_pi = %f"%(E_pi))
        
        E = E_m + E_pi
        print("E = %f"%(E))

        # change the value of the activation energy with or without a binding energy of a V-V pair, or more.
                
        if b_first_neighbours_Va >= 1 or b_second_neighbours_Va >= 1:
            rate = v*math.exp(-(E + E_bVaVa)/(kT))
        
        elif b_first_neighbours_Cu >= 1 or b_first_neighbours_Cu >= 1:
            rate = 1.0
            
        elif b_first_neighbours_Va == 0 and b_second_neighbours_Va == 0:
            rate = diffusion_rate_FeVa
            
        print ("Rate = %.1f\n"%(rate)) 
        #n = counti()
        #print("interactions = %i"%(n))
        return rate