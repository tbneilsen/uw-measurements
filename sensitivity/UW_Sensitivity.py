# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

def UW_Sensitivity(f, Xss, sens, source, mic):
    "This function calibrates the received signal \
    of a hydrophone by applying the sensitivity chart \
    of the hydrophone"
    
    hydro4034 = np.genfromtxt('4034 Hydrophone sensitivity.csv', dtype=float, delimiter=',')
    hydro4038 = np.genfromtxt('4034 Hydrophone sensitivity.csv', delimiter=',')
    project4034 = np.genfromtxt('4034 Hydrophone sensitivity.csv', delimiter=',')
    project4038 = np.genfromtxt('4034 Hydrophone sensitivity.csv', delimiter=',')
    
    if mic == 4034:
        f_mic = 1000*hydro4034[:,0]
        dB_mic = hydro4034[:,1]
        
    if mic == 4038:
        f_mic = 1000*hydro4038[:,0]
        dB_mic = hydro4038[:,1]
        
    if source == 4038:
        f_source = 1000*project4038[:,0]
        dB_source = project4038[:,1]
        
    calib = interp.interp1d(f_mic, dB_mic, kind='slinear')
    new_spec = Xss*sens*1e-9/(10**(calib(f)/20))
    
    return new_spec



plt.plot(f,test)
