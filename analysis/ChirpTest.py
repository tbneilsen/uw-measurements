# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 10:58:45 2021

@author: Corey Dobbs
"""
import numpy as np
from scipy.signal import chirp
#import TankCharacterization as TC
#import matplotlib.pyplot as plt


#times at which to evaluate the array for creating the chirp (sig)
chrp0 = 0 #chirp start time (s) (MUST START AT t=0 for CHIRP func)
chrp1 = 0.5 #chirp stop time (s)
f_0 = 1000 #Hz start freq
f_1 = 10000 #Hz end freq
fs = 48000 #sampling rate should be min = 2*f_1
trl0 = 0.5 #trailing/leading zeros
noises = True #compute a noisy signal or not
tim = np.linspace((chrp0),(chrp1),int(fs*(chrp1-chrp0)))
sig1 = chirp(tim,f_0,chrp1,f_1,method='linear')
#time array for plotting and putting in lead/trail zeros
time = np.linspace(0,(chrp1+2*trl0),int((chrp1+2*trl0)*fs))
nzeros = int(trl0*fs)
sig = np.pad(sig1,(nzeros,nzeros),'constant',constant_values=(0,0))

