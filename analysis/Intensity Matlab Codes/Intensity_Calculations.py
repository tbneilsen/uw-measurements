# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:38:29 2021

@author: Corey Dobbs
"""

#%%

#p-p method

import numpy as np
import matplotlib.pyplot as plt
import UWIntensityFuncs as uwi
import computeBlockFFT as cbFFT
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu


I_ref = 1e-12#6.61e-19
#TODO double check I_ref
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
x1 = byu.binfileload(filepath,'ID',ID_num,1)/.1
x2 = byu.binfileload(filepath,'ID',ID_num,2)/.1
    
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

#%%
import numpy as np
import matplotlib.pyplot as plt
#import UWIntensityFuncs as uwi
#import computeBlockFFT as cbFFT
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu
#particle motion detector - preliminary measurements

fs = 150000
dt = 1/fs

filepath1 = "D:/uw-acoustics-research/2021-07-15/scan1"
filepath2 = "D:/uw-acoustics-research/2021-07-15/scan2"

u1 = byu.binfileload(filepath1, 'ID', 0, 0)
u2 = byu.binfileload(filepath2, 'ID', 1, 1)

t = np.arange(300000)/fs

p1 = byu.binfileload(filepath1, 'ID', 0, 3)
p2 = byu.binfileload(filepath2, 'ID', 1, 3)

I1 = p1*u1
I1avg = np.sum(u1*p1*dt)/3

plt.figure()
plt.plot(t,u1)
plt.figure()
plt.plot(t,u2)

plt.figure()
plt.plot(t,I1)

#%% Repeat of Gabe's measurements with PMD
import numpy as np
import matplotlib.pyplot as plt
import computeBlockFFT as cbFFT
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu

fs = 150000
t = np.arange(450000)/fs
ns = 450000
df = fs/ns #sample spacing in frequency domain
fss = np.transpose(np.arange(0,(fs/2),df)) #single-sided frequency array
I_ref = 6.61e-19

#window and weighting functions
w = np.hanning(ns)
W = np.mean(w*np.conj(w)) #used to scale the ASD for energy conservation
    

#41 cm measurement
filepath = "D:/uw-acoustics-research/2021-07-16/scan8/"
ID_num = 8
#TODO cycle through ID numbers, plot pressure, and check if any show an actual signal
u_x = byu.binfileload(filepath, 'ID',0, ID_num)
u_y = byu.binfileload(filepath, 'ID',1,ID_num)
u_z = byu.binfileload(filepath, 'ID', 2, ID_num)
u_mag = np.sqrt(u_x**2 + u_y**2 + u_z**2)
p = byu.binfileload(filepath, 'ID',3, ID_num)
#TODO these values are probably in volts, we need to find the sensitivity
#p = p - np.mean(p)
#u_y = u_y - np.mean(u_y)

I_y = p*u_y
I_mag = p*u_mag

Xss, num_Blocks = cbFFT.computeBlockFFT(I_y,ns,w,W)

plt.figure(2)
plt.grid()
plt.plot(fss,10*np.log10(Xss[0]/I_ref))
plt.xlim([0,3000])
#plt.ylim([-30,110])
plt.title("PMD 41 cm from source")
plt.xlabel("Frequency (Hz)")
plt.ylabel("I (dB re 6.61e-19 W/m^2)")
