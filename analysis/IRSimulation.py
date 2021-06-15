# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 15:23:29 2021

@author: cvongsaw
"""
import ESAUdata as data
import byuarglib as byu
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab
import ESAUResponse as res
params = {'legend.fontsize': 15,
          'figure.figsize': (15, 10),
         'axes.labelsize': 24,
         'axes.titlesize':28,
         'axes.titleweight':'bold',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large',
         'lines.linewidth':2}
pylab.rcParams.update(params)


"""___SIGNALS____"""
from scipy.signal import chirp
#times at which to evaluate the array for creating the chirp (sig)
chrp0 = 0 #chirp start time (s) (MUST START AT t=0 for CHIRP func)
chrp1 = 0.5 #chirp stop time (s)
f_0 = 1e3 #Hz start freq
f_1 = 10e3 #Hz end freq
fs = 50e3 #sampling rate should be min = 2*f_1
trl0 = 0.5 #trailing/leading zeros
tim = np.linspace((chrp0),(chrp1),int(fs*(chrp1-chrp0)))
sig1 = chirp(tim,f_0,chrp1,f_1,method='linear')
#time array for plotting and putting in lead/trail zeros
time = np.linspace(0,(chrp1+2*trl0),int((chrp1+2*trl0)*fs))
nzeros = int(trl0*fs)
sig = np.pad(sig1,(nzeros,nzeros),'constant',constant_values=(0,0))
plt.figure()
plt.plot(time,sig)
plt.title(f'Sine-Swept Signal {f_0/1000}-{f_1/1000} kHz')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.grid()

####################################
"""ADD RANDOM NOISE TO THE SYSTEM"""
####################################
noise = np.random.normal(0, .1, sig.shape)
noisy = sig + noise #simple addition of noise
noisy1 = np.convolve(sig,noise,mode='same') #convolving noise into the chirp
plt.figure()
plt.plot(time,noisy)
#plt.plot(time,noisy1,'--')
plt.title('Noisy Sine-Swept Signal')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
#plt.legend(['noise added','noise convolved'])
plt.grid()


"""___Simulated Impulse Response___"""
########################################
#arbitrary impulse response for testing#
########################################
from scipy.signal import impulse, resample, unit_impulse
system = ([1.0],[1.0,2.0,1.0]) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html
timp,imp = impulse(system)#,N=int(fs*(chrp1+trl0)/2))
#make the simulated impulse response the same length as the signal by padding zeros
imp = np.pad(imp,int((len(sig)-len(imp))/2),mode='constant')
################
#delta function#
################
delta = unit_impulse(len(sig),int(0.5*fs))
#########
#Gausian#
#########
from scipy.signal.windows import gaussian
gauss = gaussian(int(len(time)),std=1)
tg = np.linspace(0,np.max(time),int(len(time)))

"""___Plot Impulsive Simulation___"""
"""___This is the target for the Response Code to Return___"""
"""___!!!___CHANGE_INPUTS_HERE___!!!___"""
#time,imp
#delta,time
#gauss,time or tg
RES = imp #IR signal to be tested
t = time #time array for the IR signal to be tested
edit = 1 #absicssa scaling error???

plt.figure()
plt.plot(t,RES)
plt.title('Arbitrary Simulated Impulse Response')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')

"""___Convolve Simulated IR w/ Chirp___"""
import scipy.signal as sci
sig_flip = np.flip(sig)
cal1 = sci.correlate(sig_flip,RES,mode='same',method='auto')
cal = np.convolve(sig,RES,mode='same')
gen = sig


"""___Obtain IR & FRF back out___"""
hsys,tsys,Hsys,fsys = res.SysResponse(cal,gen,fs,tgate=0,wiener=True,domain='f')
#tsys = tsys*1000 #convert from s to ms for plotting purposes

FRFi = np.fft.fft(RES)
Fi = np.fft.fftfreq(len(RES),d=1/(len(RES)/max(t)))
Fiss = Fi[0:int(len(Fi)/2)]/1000 #convert from Hz to kHz
FRFiss = 2*FRFi[0:(int(len(FRFi)/2))]
FRFi_dB = 10*np.log10(np.abs(FRFiss))


"""___Roll the Shape of the Time-Domain___"""
#for some reason the IR hsys is returned shifted such that the 
#alignment is a bit off so I rolled it over to align
#ro = np.gradient(hsys,1/fs)
#rol = np.where(hsys<0)
#roll = int(tsys[int(np.max(rol))]/1000*fs)
roll = int(fs*(.74999))#/1000)
hsys = np.roll(hsys,roll)
Hss = 2*Hsys[0:(int(len(Hsys)/2))]   #convert to single-sided FRF
fss = fsys[0:int(len(fsys)/2)]/1000    #convert to single-sided from Hz to kHz
Hss_dB = 10*np.log10(np.abs(Hss)) 

###############################
"""___PLOT FOR COMPARISON___"""
###############################
plt.figure() 
plt.plot(t,RES,linewidth=6)
#plt.title('Simulated Impulse Response')
plt.title('Simulated Impulse Response')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.grid()

#plt.figure()
plt.plot(tsys,np.abs(hsys),'--')
#plt.title('System Response w/ Arb. IR & Gen. Chirp')
#plt.xlabel('time (s)')
#plt.ylabel('Amplitude')
plt.legend(['Simulated IR','System Response Code'])
#plt.grid()

plt.figure()
plt.plot(Fiss*edit,FRFi_dB,linewidth=6)
#plt.title('FRF of Simulated IR from direct fft')
#plt.xlabel('Frequency (kHz)')
#plt.ylabel('Amplitude dB')
#plt.grid()

#plt.figure()
plt.plot(fss,Hss_dB,'--')
plt.title(f'Frequency Response Compar. of Chirp {f_0/1000}-{f_1/1000} kHz')
#plt.title('FRF From System Response w/ Arb. IR & Gen. Chirp')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Amplitude dB')
plt.legend(['Simulated Frequency Response','Frequency Response Code'])
plt.grid()


"""
apple = np.convolve(sig,sig_flip,mode='same')
pear = sci.correlate(sig,sig,mode='same',method='auto') 
plt.figure()
plt.plot(t,apple)
plt.plot(t,pear,'--')
plt.legend(['convolv','correlate'])
plt.title('Delta?')
plt.xlabel('time(s)')
"""