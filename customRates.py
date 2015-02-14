# here a custom rate is calculated for unique events in the simulation
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham

from KMCLib import *
import numpy
import math

#### values required for vacancy diffusion, energy in eV, T in K, v is the jump attemps in s^-1

E_m = 0.55
E_b = 0.2 #binding energy of V-V
k = 0.862e-4
T = 560
v = 1e13
kT = k*T
rate = v*math.exp(-E_m/(kT))
v_rate = 1.0


class CustomRateCalculator(KMCRateCalculatorPlugin):

    def rate(self, geometry, types_before, types_after, rate_constant, process_number, coordinate):

        # For every vacancy, check the neighbours around it to see if they are V.
        V_neighbours = len([ e for e in [types_before[1], types_before[2], types_before[3], types_before[4], types_before[5],types_before[6], types_before[7], types_before[8]] if e == "V"])
                       
        print("V_neighbours = %i \n"%(V_neighbours))
        
        #change the value of the activation energy with or without a binding energy of a V-V pair, or more.        
        if V_neighbours >= 1:
            v_rate = v*math.exp(-(E_m + E_b)/(kT))
            print("v_rate = %.1f"%(v_rate))
            
        elif V_neighbours == 0:
            v_rate = rate
        
        return v_rate

    def cutoff(self):

        #To determine the radial cutoff of the geometry around the central lattice site to cut out and send down to the custom rate function.
        #Restricts the custom rate to look in the primitive cell internal coordinates (float)

        # is this number a fuction of a?
        return 1.0