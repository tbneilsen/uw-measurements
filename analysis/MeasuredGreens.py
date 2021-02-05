# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:07:59 2020

Determining the Greens function of the tank by deconvolving or dividing out the 
generated signal response to obtain only the tank response. 

@author: cvongsaw
"""


def MeasGreen(rec,gen,fs,SysResponse,AEgirPose,RanPose,D=0.2,T=16.0,S=0.03,pressure = True):
    """
    Parameters
    ----------
    rec:    ndarray of float;
            time (unitflag = 0) or frequency (unitflag = 1) domain of the 
            received signal in pressure (if pressure = True) or Voltage (if 
            pressure = False). 
    gen:    ndarray of float;
            time (unitflag = 0) or frequency (unitflag = 1) domain of the 
            generated signal in pressure (if pressure = True) or Voltage (if 
            pressure = False). 
    fs:     float; 
            Sampling frequency
    SysResponse:    ndarray;
                    This is the frequency response of the whole measurment 
                    chain found between two close points using 
                    SystemFResponse.py SystemResponse
    AEgir_Pose: tuple; 
                AEgir TCP position (x,y,z) in the 'tank' frame 
    Ran_pose:   tuple;
                Ran TCP position (x,y,z) in the 'tank' frame
    D   : float, optional;
        water depth (m) where 0<= D <=1000m
    T   : float, optional;
        temperature in Celcius where -2<= T <=24.5 
    S   : float, optional;
        salinity where 0.030<= S <=0.042 grams salt per kg H20            
    pressure:   Boolean {True or False}; optional;
                Default (True) for providing a Pressure input and not a Voltage. 
                (False) for a Voltage input which does not take into account 
                the sensitivity of the transducers or measurement chain. 
                     
    Returns
    -------
    H_tank: ndarray of float; 
            Complex two-sided Greens function of Frequency Response of the Tank 
            envrionment 
    f:  ndarray of float;
        frequency array matching the frequency response H_tank             
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    This greens function is obviously relative to the individual points of
    Source and Receiver in the tank. To use this, you will need to also run
    the function SystemResponse to obtain the frequency response of the 
    measurement chain (transducers,etc.)
    
    last modified 9/09/2020
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from TimeGate_UnderwaterTank import timegateTank
    print('printing measurement position details...')
    print('')
    tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank')
    
    #load in recorded and generated signal at some specific AEgirPose and RanPose
    #take the fft to conver to frequency domain and obtain its frequency array
    Rec = np.fft.fft(rec)#-np.mean(rec))
    Gen = np.fft.fft(gen)+1e-12#-np.mean(gen))
    f = np.fft.fftfreq(len(rec),1/fs)#fs/(len(rec)-1))
    
    #Need to interpolate SystemResponse to match size of rec & gen lined up with 
    #the given frequency array. To do this, convert back to time domain to avoid 
    #interpolation of complex values using ifft. Pad the time domain with zeros 
    #on the end in order get back to the correct size and return the freq. 
    #domain using fft. 
    ht = np.fft.ifft(SysResponse)
    nzeros =int(np.abs((len(rec)-len(SysResponse))))
    sys = np.pad(ht,(0,nzeros),'linear_ramp',end_values=(0,0))
    Sys = np.fft.fft(sys)
    
    #Greens function or frequency response of the tank at specific AEgir & Ran
    #positions. eq. 1.7.1 Leishman 560 notes 2019.
    #Htank = Rec/(Gen*Sys) #original idea directly from notes dividing out the 
    #system response from the recorded response to have a more pure response
    #aka a deconvolution in the time domain. Adapt for potential of dividing by 
    #zeros when using complex numbers by multiplying by 1 as (conj(sys)/conj(sys))
    #and (conj(Gen)/conj(Gen)) obtaining the magnitude on bottom. 
   
    #following an adaptation of the above with Wiener Deconvolution which uses 
    #the sigma**2 value which ads a scaling factor for adjusting further for noise
    #however this is limited to only linear time-invariant signals
    lamb = 0.005 #scaling parameter
    sigma = lamb*np.max(np.abs(Gen)) #expectation or noise or SNR
    Htank = np.conj(Sys)*np.conj(Gen)*Rec/(np.abs(Gen)**2*np.abs(Sys)**2+sigma**2)
    """
    #straight forward method without the Wiener method
    Htank = np.conj(Sys)*np.conj(Gen)*Rec/np.abs(Gen)**2*np.abs(Sys)**2
    """
    #plt.figure()
    #plt.plot(np.abs(sys))
    #plt.figure()
    #plt.plot(np.abs(np.fft.ifft(Htank)))
    #plt.figure()
    #plt.plot(f,np.abs(Gen))
    #plt.figure()
    #plt.plot(f,np.abs(Rec))
    #plt.figure()
    #plt.plot(f,np.abs(Sys))
  
    return Htank,f  