# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 15:23:29 2021

@author: cvongsaw
"""
#import ESAUdata as data
#import byuarglib as byu
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
f_0 = 100e3 #Hz start freq
f_1 = 500e3 #Hz end freq
fs = 1.5e6 #sampling rate should be min = 2*f_1
trl0 = 0.1 #trailing & leading zeros
noises = True #compute a noisy signal or not
tim = np.linspace((chrp0),(chrp1),int(fs*(chrp1-chrp0)))
sig1 = chirp(tim,f_0,chrp1,f_1,method='linear')
#time array for plotting and putting in lead/trail zeros
time = np.linspace(0,(chrp1+2*trl0),int((chrp1+2*trl0)*fs))
nzeros = int(trl0*fs)
sig = np.pad(sig1,(nzeros,nzeros),'constant',constant_values=(0,0))

#divide by convert to change Hz to kHz if 1000 or leave as Hz if 1
if f_0>=1e3:
    convert = 1000
if f_1<=1e3:
    convert = 1 

plt.figure()
plt.plot(time,sig)
#if convert == 1000:
    #plt.title(f'Swept-Sine Signal {f_0/1000}-{f_1/1000} kHz')
#if convert == 1:
    #plt.title(f'Swept-Sine Signal {f_0}-{f_1} Hz')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.savefig('Figure1')
plt.close()

"""ADD RANDOM NOISE TO THE SYSTEM"""
if noises == True:
    noise = np.random.normal(0, .1, sig.shape)
    noisy = sig + noise #simple addition of noise
    noisy1 = np.convolve(sig,noise,mode='same') #convolve noise into the chirp
    plt.figure()
    plt.plot(time,noisy)
    #if convert == 1000:
    #    plt.title(f'Noisy Swept-Sine Signal {f_0/1000}-{f_1/1000} kHz')
    #if convert == 1:
    #    plt.title(f'Noisy Swept-Sine Signal {f_0}-{f_1} Hz')
    plt.xlabel('time (s)')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.savefig('Figure2')
    plt.close()

"""___Simulated Impulse Response___"""
########################################
#arbitrary impulse response for testing#
########################################
from scipy.signal import impulse, resample, unit_impulse
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html
system = ([1.0],[3.0,2.0,1.0]) 
timp,imp = impulse(system)

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


"""Impulse Simulation with REFLECTION"""
#append 0.5*imp to imp to simulate a reflection? to be timegated 




"""___Plot Impulsive Simulation___"""
"""___This is the target for the Response Code to Return___"""
"""___!!!___CHANGE_INPUTS_HERE_When_Other_IR_Desired___!!!___"""
#time,imp
#delta,time
#gauss,time or tg
RES = imp #IR signal to be tested
t = time #time array for the IR signal to be tested
edit = 1 #absicssa scaling error???

plt.figure()
plt.plot(t,np.abs(RES))
plt.title('Simulated Impulse Response (IR)')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.savefig('Figure3')
plt.close()






"""___Convolve Simulated IR w/ Chirp___"""
#import scipy.signal as sci
#sig_flip = np.flip(sig)
#cal1 = sci.correlate(sig_flip,RES,mode='same',method='auto')
cal = np.convolve(sig,RES,mode='same')
gen = sig
plt.figure()
plt.plot(t,cal)
#if convert == 1000:
#    plt.title(f'{f_0/1000}-{f_1/1000} kHz Chirp Convolved w/ Simulated IR')
#if convert == 1:
#    plt.title(f'{f_0}-{f_1} Hz Chirp Convolved w/ Simulated IR')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.savefig('Figure4')
plt.close()

if noises == True:
    calnoise = np.convolve(noisy,RES,mode='same')
    plt.figure()
    plt.plot(t,calnoise)
    #if convert == 1000:
    #    plt.title(f'{f_0/1000}-{f_1/1000} kHz Noisy Chirp Convolved w/ Simulated IR')
    #if convert == 1:
    #    plt.title(f'{f_0}-{f_1} Hz Noisy Chirp Convolved w/ Simulated IR')
    plt.xlabel('time (s)')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.savefig('Figure5')
    plt.close()

"""___Obtain IR & FRF back out___"""
hsys,tsys,Hsys,fsys = res.SysResponse(cal,gen,fs,tgate=0,wiener=True,domain='f')

FRFi = np.fft.fft(RES)
Fi = np.fft.fftfreq(len(RES),d=1/(len(RES)/max(t)))
Fiss = Fi[0:int(len(Fi)/2)]/convert #convert from Hz to kHz
FRFiss = 2*FRFi[0:(int(len(FRFi)/2))]
FRFi_dB = 10*np.log10(np.abs(FRFiss))

"""___Roll the Shape of the Time-Domain___"""
#The IR is shifted to the end of the array, such that the tail
#spills over to the beginning of the ray. The array must be rolled
#for alignment w/ the actual IR. However, the number of zeros must
#be equal for both leading and trailing zeros. 
roll = int(0.5*len(hsys)-1)
hsys = np.roll(hsys,roll)
Hss = 2*Hsys[0:(int(len(Hsys)/2))]   #convert to single-sided FRF
fss = fsys[0:int(len(fsys)/2)]/convert    #convert to single-sided from Hz to kHz
Hss_dB = 10*np.log10(np.abs(Hss)) 


if noises == True:
    #NOISY VERSION OF SYSTEM RESPONSE
    hsysn,tsysn,Hsysn,fsysn = res.SysResponse(cal,gen,fs,tgate=0,wiener=True,domain='f')
    """___Roll the Shape of the Time-Domain___"""
    #The IR is shifted to the end of the array, such that the tail
    #spills over to the beginning of the ray. The array must be rolled
    #for alignment w/ the actual IR. However, the number of zeros must
    #be equal for both leading and trailing zeros. 
    hsysn = np.roll(hsysn,roll)
    Hssn = 2*Hsysn[0:(int(len(Hsysn)/2))]   #convert to single-sided FRF
    fssn = fsysn[0:int(len(fsysn)/2)]/convert    #convert to single-sided from Hz to kHz
    Hss_dBn = 10*np.log10(np.abs(Hssn)) 



"""___PLOT FOR COMPARISON___"""
plt.figure() 
plt.plot(t,np.abs(RES),linewidth=6)
plt.plot(tsys,np.abs(hsys),'--',linewidth=3)
#if convert == 1000:
#    plt.title(f'{f_0/1000}-{f_1/1000} kHz Impulse Response w/ fs={fs/1000}kHz')
#if convert == 1:
#    plt.title(f'{f_0}-{f_1} Hz Impulse Response w/ fs={fs}Hz')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')
plt.legend(['Simulated IR','Deconvolved IR'])
plt.grid()
plt.savefig('Figure6')
plt.close()

if noises == True:
    plt.figure() 
    plt.plot(t,np.abs(RES),linewidth=6)
    plt.plot(tsysn,np.abs(hsysn),'--',linewidth=3)
    #if convert == 1000:
    #    plt.title(f'{f_0/1000}-{f_1/1000} kHz Impulse Response w/ Noise & fs={fs/1000}kHz')
    #if convert == 1:
    #    plt.title(f'{f_0}-{f_1} Hz Impulse Response w/ Noise & fs={fs}Hz')
    plt.xlabel('time (s)')
    plt.ylabel('Amplitude')
    plt.legend(['Simulated IR','Deconvolved IR'])
    plt.grid()
    plt.savefig('Figure7')
    plt.close()

plt.figure()
plt.plot(Fiss*edit,FRFi_dB,linewidth=6)
plt.plot(fss,Hss_dB,'--',linewidth=3)
#if convert == 1000:
#    plt.title(f'Frequency Response of {f_0/1000}-{f_1/1000} kHz Signal')
#    plt.xlabel('Frequency (kHz)')
#if convert == 1:
#    plt.title(f'Frequency Response of {f_0}-{f_1} Hz Signal')
#    plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude dB')
plt.legend(['Simulated Frequency Response','Deconvolved Frequency Response'])
plt.grid()
buffer_limit = f_1+(f_1-f_0)*0.01
#plt.xlim(f_0-buffer_limit,f_1+buffer_limit)
plt.savefig('Figure8')
plt.close()

if noises == True:
    plt.figure()
    plt.plot(Fiss*edit,FRFi_dB,linewidth=6)
    plt.plot(fssn,Hss_dBn,'--',linewidth=3)
    #if convert == 1000:
    #    plt.title(f'Frequency Response of {f_0/1000}-{f_1/1000} kHz Noisy Signal')
    #    plt.xlabel('Frequency (kHz)')
    #if convert == 1:
    #    plt.title(f'Frequency Response of {f_0}-{f_1} Hz Noisy Signal')
    #    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude dB')
    plt.legend(['Simulated Frequency Response','Deconvolved Frequency Response'])
    plt.grid()
    buffer_limit = f_1+(f_1-f_0)*0.01
    #plt.xlim(f_0-buffer_limit,f_1+buffer_limit)
    plt.savefig('Figure9')
    plt.close()