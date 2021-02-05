# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:33:19 2020

@author: cvongsaw
"""

"""
#####################################################
#this is test information that can be removed later.# 
import numpy as np                                  #
fs = 1e6                                            #
A = (0.5,3,0.15)                                    #
R = (0.5,1,0.15)                                    #
t = np.arange(0,1,1/fs)                             #
rec = np.zeros(len(t))                              #
gen = np.zeros(len(t))                              #
rec[50:59]=1    #delayed impulse mimic received     #
gen[0:9]=1      #generated impulse                  #
#####################################################
"""

def SystemResponse(cal,gen,fs,AEgirPose,RanPose,tstart=0.4499,D=0.2,T=16.0,S=0.03,pressure = True):
    """
    Parameters
    ----------
    cal:    ndarray of float;
            time (unitflag = 0) or frequency (unitflag = 1) domain of the 
            received signal in pressure (if pressure = True) or Voltage (if 
            pressure = False). 
    gen:    ndarray of float;
            time (unitflag = 0) or frequency (unitflag = 1) domain of the 
            generated signal in pressure (if pressure = True) or Voltage (if 
            pressure = False). 
    fs:     float; 
            Sampling frequency
    AEgir_Pose: tuple; 
                AEgir TCP position (x,y,z) in the 'tank' frame 
    Ran_pose:   tuple;
                Ran TCP position (x,y,z) in the 'tank' frame
    tstart:     float;
                leading zeros time delay that should be expected in the signal. 
                this is a setting from ESAU when generating a signal. 
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
    Hf:     ndarray of float 
            Complex two-sided Frequency Response to account for all transducer,
            amplifier, etc. in the measurement chain.     
    f:  ndarray of float
        frequency array
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Measurements should be taken with source and receiver close together and 
    in the center of the tank so it is easy to time gate the signal. Not too 
    close that nonlinear effects occur. The callibration position should already
    be hard coded into ESAU. 
    
    The calibration measurement should be recorded on ch0. This code is currently
    only setup to obtain the Through The System Response relative to the gen
    signal and the calibration measurement input into ESAU via the ch0 card input.
    
    last modified 11/16/2020
    """
    #TimeGate_UnderwaterTank.py timegateTank
    from TimeGate_UnderwaterTank import timegateTank
    import numpy as np
    import matplotlib.pyplot as plt
    
    print('')
    print('printing calibration position details...')
    print('')
    tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank') 
    
    #convert leading zeros time length to samples of leading zeros
    start = int(tstart*fs)
    fin = int((tgate+tstart)*fs)
    #gate the signal for only samples from start to finish 
    gen_gate = gen[start:fin]
    cal_gate = cal[start:fin]
    #plt.figure()
    #plt.plot(gen_gate)
    #plt.title('gen_gate')
    #plt.figure()
    #plt.plot(cal_gate)
    #plt.title('cal_gate')
    
    #convert gated signal to freq domain
    Cal = np.fft.fft(cal_gate)#-np.mean(cal_gate))
    Gen = np.fft.fft(gen_gate)#-np.mean(gen_gate))
    
    #frequency array using speed(c) from TimeGate_UnderwaterTank.py calculation
    f = np.fft.fftfreq(len(Cal),d=1/fs)
    k = 2*np.pi*f/c
    #remove phase shift due to propagation through distance(dis) in water
    Cal_mod = Cal[:,0]*np.exp(-1j*k*dis)
    
    #Determine the Frequency Response of the System (transducers,amps, etc.)
    #Leishman notes Phscs560 eq 1.7.2 while also multiplying by 1 aka 
    #conj(gen)/conj(gen) in order to remove any zeros from complex values used.
    #When converting back to the time domain with an ifft to obtain an impulse 
    #response, this looks similar to a Cross correlation operation according to
    #Brian Anderson but not exactly. 
    
    #following an adaptation of the above with Wiener Deconvolution which uses 
    #the sigma**2 value which ads a scaling factor for adjusting further for noise
    #however this is limited to only linear time-invariant signals
    ###not sure this should be done for the system response
    lamb = 0.005 #scaling parameter
    sigma = lamb*np.max(np.abs(Gen)) #expectation or noise or SNR
    Hf = np.conj(Gen)*Cal_mod/(np.abs(Gen)**2+sigma**2)
    ht = np.fft.ifft(Hf)
    """
    #direct method of deconvolution division in freq domain without Weiner adjustment
    Hf = np.conj(Gen)*Cal_mod/(np.abs(Gen)**2)
    ht = np.fft.ifft(Hf)
    """
    #plt.figure()
    #plt.plot(Gen)
    #plt.figure()
    #plt.plot(Hf)
    #plt.title('Sys FRF')
    plt.figure()
    plt.plot(np.abs(ht))
    plt.title('sys IR')
    return Hf,f