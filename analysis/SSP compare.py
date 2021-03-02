# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 22:08:22 2020

This code was written to compare three well known underwater sound speed 
formulae as cited below. I found that they were similar and that Garrett was 
well known, accepted, and reasonable enough for our models. There are more 
computationally heavy models used however, but it was not worth my time to 
get into those in depth. 


@author: Cameron Vongsawad
"""

import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0,50,500)  #temperature in celcius
s = np.linspace(0,45,500)  #salinity grams of salt per kg of water. aka parts per 1000
d = np.linspace(0,600,500) #depth in meters


#garrett
C1 = lambda T,S,D: 1493 + 3*(T - 10) - 0.006*(T - 10)**2 - 0.04*(T - 18)**2 + 1.2*(S - 35) - 0.01*(T - 18)*(S - 35) + D/61  
#medwin & Kuperman Encyclo. of Ocean Sciences 2nd ed. 2001
C2 = lambda T,S,D: 1449.2 + 4.6*T - 0.055*T**2 + 0.00029*T**3 + (1.34 - 0.010*T)*(S - 35) + 0.016*D
#Christ, WenliSr., The ROV Manual 2nd ed. 2014 from simplified Wilson's 1960
#S is in PSU which is basically equivalent to ppt
C3 = lambda T,S,D: 1449 + 4.6*T - 0.055*T**2 + 0.0003*T**3 + 1.39*(S - 35) + 0.017*D

c1 = np.zeros(500)
c2 = np.zeros(500)
c3 = np.zeros(500)

p = 0
#D = 0.5
s = np.zeros(500)
for T in t:
    for S in s:
        for D in d: 
            c1[p] = C1(T,S,D)
            c2[p] = C2(T,S,D)
            c3[p] = C3(T,S,D)
            if p == 500:
                p = p -1
            else: 
                p = p +1

            
plt.figure()
plt.plot(c1)
plt.plot(c2)
plt.plot(c3)
plt.legend(['Garrett','Medwin & Kuperman','Wilson'])
plt.show()