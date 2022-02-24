# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:38:29 2021

@author: Corey Dobbs
"""

#%%

#p-p method
#This is a transcription of Gabe's thesis code. 

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

u_x = byu.binfileload(filepath, 'ID',0, ID_num)
u_y = byu.binfileload(filepath, 'ID',1,ID_num)
u_z = byu.binfileload(filepath, 'ID', 2, ID_num)
u_mag = np.sqrt(u_x**2 + u_y**2 + u_z**2)
p = byu.binfileload(filepath, 'ID',3, ID_num)
I_y = p*u_y
I_mag = p*u_mag

Xss, num_Blocks = cbFFT.computeBlockFFT(I_y,ns,w,W)

plt.figure(2)
plt.grid()
#plt.plot(fss,10*np.log10(Xss[0]/I_ref))
plt.plot(t,p)
#plt.xlim([0,3000])
#plt.ylim([-30,110])
plt.title("PMD 41 cm from source")
#plt.xlabel("Frequency (Hz)")
#plt.ylabel("I (dB re 6.61e-19 W/m^2)")


#%%
#Test OMNI sensor
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu

import sys
sys.path.append("D:/uw-acoustics-research/personal-codes/analysis/")
import UWIntensityFuncs as uwi
plt.close('all')


fs = 150000
t = np.arange(450000)/fs

filepath = "D:/uw-acoustics-research/2021-09-08/"

rho = 1026
c = 1471
sig_cond_bk = .01 #V/Pa

ux = byu.binfileload(filepath, 'ID', 2, 0)
uy = byu.binfileload(filepath, 'ID',2,1)
uz = byu.binfileload(filepath,'ID',2,2)
p_pmd1 = byu.binfileload(filepath,'ID',2,3) #OMNI sensor plugged in with velocity channels

p_bk = byu.binfileload(filepath,'ID',3,1)/sig_cond_bk
p_pmd=byu.binfileload(filepath,'ID',3,2)

ux1 = byu.binfileload(filepath,'ID',1,0)
uy1 = byu.binfileload(filepath,'ID',1,1)
uz1 = byu.binfileload(filepath,'ID',1,2)


# replacing zero values
ux = np.where(ux>0,ux,1e-16)
uy = np.where(uy>0,uy,1e-16)
uz = np.where(uz>0,uz,1e-16)
p_pmd = np.where(p_pmd>0,p_pmd,1e-16)


vx,vy,vz,omni = uwi.PMD_processing(fs,ux,uy,uz,p_pmd)



#Find frequency spectrum of OMNI sensor

#TODO check monitor signals between these and Gabe's thesis signals
#Compare OASPL, check frequency range, ESAU settings, conditioner box settings(10mV/Pa?) 
pref=1e-6
p_pmd = 1e-6*10**((20*np.log10(abs(p_pmd)))/20)

spec,fss,OASPL = byu.autospec(p_pmd-np.mean(p_pmd),fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
spec1,fss1,OASPL1 = byu.autospec(p_bk-np.mean(p_bk),fs,ns=2**15,N=0,unitflag=0,pref=1e-6)

plt.figure()
plt.plot(fss,20*np.log10(spec/pref))
plt.grid()
plt.title('Freq Spectrum of OMNI sensor')
plt.xlim([0, 10000])

#Find frequency spectra of x-direction sensors
specx,fssx,OASPLx = byu.autospec(ux,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
specy,fssy,OASPLy = byu.autospec(uy,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
specz,fssz,OASPLz = byu.autospec(uz,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
plt.figure()
plt.semilogx(fssx,20*np.log10(np.abs(specx)/pref),fssy,20*np.log10(np.abs(specy)/pref),fssz,20*np.log10(np.abs(specz)/pref))
plt.grid()
plt.title('Freq Spectrum of v-Sensor')
plt.legend(['x','y','z'])
#plt.xlim([0, 10000])
#plt.ylim([0,0.25e-14])

#Plotting
plt.figure()
plt.plot(fss,vx,fss,vy,fss,vz)
plt.grid()
plt.legend(('x','y','z','p'))
# plt.xlim([10, 10000])
# plt.ylim([0,1e-7])

plt.figure()
plt.plot(fss,omni)
plt.xlim([0,10000])
plt.ylim([0,100])
plt.grid()
plt.xlabel('Frequency (Hz)')
plt.ylabel('Pressure (Pa)')

#Attempt at calculating time-averaged intensity
T = ns/fs
dt = 1/fs
v_mag = np.sqrt(ux**2 + uy**2 + uz**2)
Imag = np.sum(p_bk*v_mag*dt)/T

#%% code to generate white noise
import generator 
plt.close('all')
fs = 150000
N = 2*fs
x = generator.white(N, state=None)
trl0 = 0.5
nzeros = int(trl0*fs)
sig = np.pad(x,(nzeros,nzeros),'constant',constant_values=(0,0))

from scipy.signal import butter,filtfilt

def butter_lowpass_filter(data, cutoff, fs, order):
    nyq = 0.5*fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

sig = butter_lowpass_filter(sig, 10000, fs, 10)

plt.figure()
plt.plot(sig)

spec,fss,OASPL = byu.autospec(x,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
plt.figure()
plt.plot(fss,spec)

#np.savetxt("WhiteNoise1.txt",sig)







