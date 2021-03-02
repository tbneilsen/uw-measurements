# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:39:55 2021

@author: cvongsaw
"""

def IR(rec,gen,fs):
    """
    Parameters
    ----------
    rec:    ndarray of float of size 1;
            timedomain of the received signal in pressure (if pressure = True) 
            or Voltage (if pressure = False). 
    gen:    ndarray of float;
            time domain of the generated signal in pressure (if pressure = True) 
            or Voltage (if pressure = False). 
    fs:     float; 
            Sampling frequency in Hz
       
    Returns
    -------
    ht:     ndarray of float;
            Real impulse response (IR) of the measurement chain neglecting 
            effects of the water and tank environment through timegating only 
            direct signal with small propagation assumption. 
    
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    The IR h(t) is determined following Gemba(2014) eq. 3.3.6 scaling similar 
    to a matched filter and deconvolving in order to obtain the pure delay of h(t).
    Cite: "Characterization of underwater acoustic sources recorded in reverberant
    environments with application to scuba signatures" Gemba (2014) Dissertation
        
    last modified 02/16/2021
    """
    import numpy as np
    import scipy.signal as sci 
        
    #Gemba (2014) dissertation eq 3.3.3, 3.3.5 and 3.3.6
    #Get the temporal inverse of the generated signal so the xcorr is 
    #basically a convolution in the time domain. The temporal inversion also 
    #makes sure there is a phase inversion and thus the xcorr results in a 
    #pure delay of h(t)
    #rec(t)*conj(gen(-t)) = h(t)*gen(t)*conj(gen(-t)) eq 3.3.5 solve for h(t)
    #this method of xcorr is equivalent to a Convolution since a convolution 
    #takes the temperal inverse of the array (inverse filter).
    gen_flip = np.flip(gen)
    #xcorr of the recorded and temporal inverse generated (xcorr takes conj of
    #the gen_flip as part of its process) acorr is autocorrelation of the gen
    #such that xcorr = rec(t)*conj(gen(-t)) from above equation and 
    #acorr = gen(t)*conj(gen(-t)) 
    acorr = sci.correlate(gen,gen_flip,mode='full',method='auto')
    xcorr = sci.correlate(rec,gen_flip,mode='full',method='auto')
    #It is necessary to take these two parts into the frequency domain in order
    #to divide them out to solve for h(t). This division in the frequency domain
    #is equivalent to a deconvolution in the time domain. which gives H(f)
    #and then the ifft(H(f)) = h(t)
    #Frequency Domain
    Xcorr = np.fft.fft(xcorr)
    Acorr = np.fft.fft(acorr)
    #Single Sided
    Xcorr = 2*Xcorr[0:(int(len(Xcorr)/2))] #******edit normalization
    Acorr = 2*Acorr[0:(int(len(Acorr)/2))]
    
    #Determine Frequency band of the gen signal to apply a bandpass filter to 
    #Acorr so as to not divide by various close to zero values causing significant
    #high frequency processing noise. This cuts the leading and trailing zeros out
    #to not obtain noise outside of the signal when performing the Deconvolution. 
    # TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO 
    cond = 0.0001*np.abs(np.max(Xcorr))
    band = np.argwhere(np.abs(Xcorr)>cond)
    hpass = min(band)
    lpass = max(band)
    #deconvolution interval
    Xd = Xcorr[int(hpass):int(lpass)]
    Ad = Acorr[int(hpass):int(lpass)]
    Deconv = Xd/Ad
    
    #Wiener Deconvolution deals with the near zero values which cause processing noise
    #using WDeconv until I figure out how to determine start/stop idx of bandpass for 
    #any arbitrary signal sent into this function. 
    lamb = 0.005 #semi arbitrary? scaling parameter
    sigma = lamb*np.max(np.abs(Ad)) #expectation or noise or SNR
    WDeconv = np.conj(Ad)*Xd/(np.abs(Ad)**2+sigma**2)
    #bring back to time domain with inverse fast fourier transform
    ht = np.fft.ifft(WDeconv)
    
    return ht
    
    






    
def SysResponse(cal,gen,fs,AEgirPose,RanPose,tstart=0.5,D=0.2,T=16.0,S=0.03,pressure = True):
    """
    Parameters
    ----------
    cal:    ndarray of float;
            time domain of the received signal in pressure (if pressure = True) 
            or Voltage (if pressure = False). 
    gen:    ndarray of float;
            time domain of the generated signal in pressure (if pressure = True) 
            or Voltage (if pressure = False). 
    fs:     float; 
            Sampling frequency
    AEgir_Pose: tuple; 
                AEgir TCP position (x,y,z) in the 'tank' frame 
    Ran_pose:   tuple;
                Ran TCP position (x,y,z) in the 'tank' frame
    tstart:     float;
                leading zeros time delay that should be expected in the signal. 
                this is a setting from ESAU when generating a signal. Defaults to 
                0.5s which is a common number of leading zeros when this is written
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
    ht:     ndarray of float;
            Real impulse response (IR) of the measurement chain neglecting 
            effects of the water and tank environment through timegating only 
            direct signal with small propagation assumption. 
    t:      ndarray of float;
            Time array for the IR h(t)
    Hf:     ndarray of complex; 
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
    be hard coded into ESAU but can be changed manually. 
    
    The calibration measurement should be recorded on ch0. This code is currently
    only setup to obtain the Through The System Response relative to the gen
    signal and the calibration measurement input into ESAU via the ch0 card input.
    
    The IR h(t) is determined first by timegating the signal for only direct sound
    then following Gemba(2014) eq. 3.3.6 scaling similar to a matched filter
    and deconvolving in order to obtain the pure delay of h(t).
    Cite: "Characterization of underwater acoustic sources recorded in reverberant
    environments with application to scuba signatures" Gemba (2014) Dissertation
        
    last modified 02/16/2021
    """
    
    from TimeGate_UnderwaterTank import timegateTank
    import numpy as np
    from ESAUResponse import IR
    
    tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank') 
    ht = IR(cal,gen,fs)
    t = np.linspace(0,len(ht)/fs,len(ht))
    #convert time in sec to time in ms
    t = t*1000
    
    #convert leading zeros time length to samples of leading zeros
    #take only majority % of the zeros allowing to view just before the signal 
    #rises with % assumed roughly by fs values for proper allowance of signal.
    #but realistically these values are mostly arbitrarily trying to approach
    #very close to the initial signal considering fs.0.99999 may be enough?????
    if fs <= 5e5:
        percent = 0.99999
    if fs > 5e5 and fs <= 1e6:
        percent = 0.999
    if fs > 1e6:
        percent = 0.99
    
    #################################################################
    #********????????????????????? TO DO ***********************
    #TIMEGATE ht SIGNAL HERE!!! before obtaining Hf
    tstart = tstart *percent 
    start = int(tstart*fs)
    fin = int((tgate+tstart)*fs)
    start = 0
    fin = int(tgate*fs)
    ht = ht[start:fin]
    t = t[start:fin]
    
    f = np.fft.fftfreq(len(ht),d=1/fs)
    Hf = np.fft.fft(ht)
    Hss = 2*Hf[0:(int(len(Hf)/2))]
    fss = f[0:int(len(f)/2)]
    

    return ht,t,Hss,fss




def TankResponse(rec,gen,fs,SysFRF,pressure = True):
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
    SysFRF: ndarray;
            This is the frequency response of the whole measurment 
            chain found between two close points using 
            SystemFResponse.py SystemResponse        
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
    from ESAUResponse import IR 
    
    #obtain the IR of the recorded signal relative to the generated signal
    ht = IR(rec,gen,fs)
    #convert to Frequency domain in order to later deconvolve the system Response
    Hf = np.fft.fft(ht)
    f = np.fft.fftfreq(len(ht),d=1/fs)
    
    #Need to interpolate SystemResponse to match size of rec & gen lined up with 
    #the given frequency array. To do this, convert back to time domain to avoid 
    #interpolation of complex values using ifft. Pad the time domain with zeros 
    #on the end in order get back to the correct size and return the freq. 
    #domain using fft. 
    ht = np.fft.ifft(SysFRF)
    nzeros =int(np.abs((len(rec)-len(SysFRF))))
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
    sigma = lamb*np.max(np.abs(Sys)) #expectation or noise or SNR
    Htank = np.conj(Sys)*Hf/(np.abs(Sys)**2+sigma**2)

    return Htank,f  





