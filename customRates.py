# here a custom rate is calculated for unique events in the simulation

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

    def rate(self, geometry, elements_before, elements_after, rate_constant, process_number, coordinate):
        

        # For every vacancy, check the neighbours around it to see if they are V.
        V_neighbours = len([ e for e in [elements_before[1], elements_before[2], elements_before[3], elements_before[4], elements_before[5],elements_before[6], elements_before[7], elements_before[8]] if e == "V"])
        
        #Fe_neighbours = len([ e for e in [elements_before[1], elements_before[2], elements_before[3], elements_before[4], elements_before[5],elements_before[6], elements_before[7], elements_before[8]] if e == "Fe"])
        
        #V_neighbours = 1
       
        #print("V_neighbours = %i \n"%(V_neighbours))
        #print("Fe_neighbours = %f \n"%(Fe_neighbours))
        
        #change the value of the activation energy with or without a binding energy of a V-V pair, or more.
        
        if V_neighbours >= 1:
            v_rate = v*math.exp(-(E_m + E_b)/(kT))
            #print("v_rate = %.1f"%(v_rate))
            
        elif V_neighbours == 0:
            v_rate = rate
        
        #print("v_rate = %f"%(v_rate))
        return v_rate

    def cutoff(self):
        return 1.0