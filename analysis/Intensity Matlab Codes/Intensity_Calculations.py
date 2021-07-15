# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:38:29 2021

@author: Corey Dobbs
"""
import numpy as np
import matplotlib.pyplot as plt
import UWIntensityFuncs as uwi
import computeBlockFFT as cbFFT
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu


I_ref = 6.61e-19
L = 3
fs = 150000
ns =  L*fs #number of samples
df = fs/ns #sample spacing in frequency domain
fss = np.transpose(np.arange(0,(fs/2),df)) #single-sided frequency array
    
#window and weighting functions
w = np.hanning(ns)
W = np.mean(w*np.conj(w)) #used to scale the ASD for energy conservation
    
# Load in data, make sure binfileload.py is downloaded from byuarg library.
#This is also a good place to account for the amplification factor from 
#the NEXUS conditioning amplifier (in this case divide by .10 to account 
#for the 10 mV/Pa amplification).

filepath = "D:/uw-acoustics-research/2021-02-17_scan8"
ID_num = 2
x1 = byu.binfileload(filepath,'ID',ID_num,1)/.10
x2 = byu.binfileload(filepath,'ID',ID_num,2)/.10
    
Xss1, num_Blocks = cbFFT.computeBlockFFT(x1,ns,w,W)
Xss2, num_Blocks = cbFFT.computeBlockFFT(x2,ns,w,W)
    
    
Xss = np.concatenate((np.reshape(Xss1,(1,np.size(Xss1,0),np.size(Xss1,1))),np.reshape(Xss2,(1,np.size(Xss1,0),np.size(Xss1,1))))) 
    
probe_config = np.array([[0,.17,0],[0,-.17,0]])
        
I, Q, U, T, E, EL, I_mag, I_dir, Q_mag, Q_dir, p0, grad_p, u, z = uwi.TRAD_func(fss,Xss,probe_config)

I_trad_y = 10*np.log10(abs(I[1,:])/I_ref)
    
Q_trad_y = 10*np.log10(abs(Q[1,:])/I_ref)

trad_mag = 10*np.log10(np.sqrt(Q[1,:]**2 + I[1,:]**2)/I_ref)

plt.figure(1)
plt.grid()
plt.plot(fss,I_trad_y,color = 'green')
plt.plot(fss,Q_trad_y,color = 'red')
plt.plot(fss,trad_mag, color = 'blue')
plt.legend(["Active Intensity","Reactive Intensity","Absolute Magnitude"])
plt.title("Y-mount 41 cm From Source")
plt.xlim([0,3000])
plt.xlabel('Frequency (Hz)')
plt.ylabel('I (dB re. 1 \u03BCPa)')
