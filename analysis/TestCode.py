# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:39:55 2021

@author: cvongsaw
"""

def IR(rec,gen,fs,wiener=False,domain='f'):
    """
    Parameters
    ----------
    rec:        ndarray of float of size 1;
                time domain of the received signal. Should be real valued.
    gen:        ndarray of float;
                time domain of the generated signal. Should be real valued.
    fs:         float; 
                Sampling frequency in Hz
    wiener:     Boolean {True or False}; optional;
                False (default) for using direct deconvolution instead of Wiener
                deconvolution in frequency domain. If (True), the Wiener 
                deconvolution is performed. Wiener deconvolution acts as a 
                regularization which helps prevent dividing by zero allowing for 
                a more robust deconvolution while maintaining an account for any 
                system response. 
    domain:     string, Optional;
                Choice of domain performs the inverse filter in the initial step
                in either the temporal domain ('t' or 'time' or 'temporal') or 
                in the frequency domain (default) ('f' or 'freq' or 'frequency') 
                which is equivalent to determining the the cross-spectral density 
                and the auto-spectral density for the use in the deconvolution. 
                The end deconvolution always occurs by division in frequency domain.
                
    Returns
    -------
    ht:         ndarray of float;
                Real valued impulse response (IR) of a measurement.
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    The IR h(t) is determined following Gemba(2014) eq. 3.3.6 scaling similar 
    to a matched filter and deconvolving in order to obtain the pure delay of h(t).
    Cite: "Characterization of underwater acoustic sources recorded in reverberant
    environments with application to scuba signatures" Gemba (2014) Dissertation.
    Also see eq. 3.3.3 and 3.3.5
    
    This also follows the directions from Farina 2000 and Farina 2007 on IR
    from swept-sines. 
    
    Also see eq. 1.7.1, 1.7.2, and 1.7.3 from Leishman 560 notes 2019.
    
    Deconvolution all in the frequency domain should be much faster computationally.
    
    Dr. Brian Anderson published a paper discussing Wiener deconvolution as a
    regularization parameter for deconvolution. He particularly discusses 
    optimizing lambda."Time reversal focusing of high amplitude sound in a 
    reverberation chamber" (2018) Willardson, Anderson, Young, Denison, Patchett.
    https://doi.org/10.1121/1.5023351
    
    last modified 5/17/2021
    """
    import numpy as np
    if domain == 'time' or 't' or 'temporal':    
        #rec(t)*gen(-t)) = h(t)*gen(t)*gen(-t) eq 3.3.5 solve for h(t)
        gen_flip = np.flip(gen) 
        #np.convolve(gen,gen_flip) == sci.correlate(gen,gen) by def.
        #the inverse filter of the function is np.convolve(gen,gen_flip)
        acorr = np.convolve(gen,gen_flip,mode='same')
        xcorr = np.convolve(rec,gen_flip,mode='same')
        #import scipy.signal as sci 
        #acorr = sci.correlate(gen,gen,mode='same',method='auto')
        #xcorr = sci.correlate(rec,gen,mode='same',method='auto')
        
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(np.abs(acorr))
        plt.title('Delta Function as result of the Inverse Filter Convolution')
        plt.xlabel('time (Samples)')
        plt.ylabel('Amplitude')
        plt.grid()
        
        #Division in the frequency domain is a deconvolution in the time domain. 
        #which gives H(f) and then ifft(H(f))=h(t)
        Xcorr = np.fft.fft(xcorr) #double-sided frequency response
        Acorr = np.fft.fft(acorr) #double-sided frequency response
        #f = np.fft.fftfreq(len(xcorr),d=1/fs)
        
    if domain == 'frequency' or 'freq' or 'f':
        #COMPUTE ALL of the deconvolution in FREQ DOMAIN instead of time domain
        Xcorr = np.fft.fft(rec)*np.conj(np.fft.fft(gen))
        Acorr = np.fft.fft(gen)*np.conj(np.fft.fft(gen))
        #f = np.fft.fftfreq(len(gen),d=1/fs)
    
    if wiener == True:
        print('Performing deconvolution via Wiener deconvolution preventing dividing by zero')
        #Wiener Deconvolution deals with the near zero values which cause processing noise
        #and high frequency aliasing.
        lamb = 0.005 #scaling parameter arbitrarily chosen 
        sigma = lamb*np.mean(np.abs(Acorr)) #expectation or noise or SNR
        WDeconv = np.conj(Acorr)*Xcorr/(np.abs(Acorr)**2+sigma**2)
        Deconv = WDeconv
    else: 
        print('Performing deconvolution via direct division in frequency domain')
        #Perform standard deconvolution by direct division in frequency domain. 
        Deconv = Xcorr/Acorr
        
    #bring back to time domain with inverse fast fourier transform (IFFT)
    ht = np.real(np.fft.ifft(Deconv)) #ensure real valued as it should be
    return ht
    
    





    
def SysResponse(cal,gen,fs,tgate=0,wiener=False,domain='f'):
    """
    Parameters
    ----------
    cal:        ndarray of float;
                Received calibration signal. Should be real valued.
    gen:        ndarray of float;
                Pure generated signal. Should be real valued.
    fs:         float; 
                Sampling frequency (Hz).
    tgate:      float,Optional;
                Time of the first wall reflection determined through timegateTank. 
                This is the time we will use to determine the time of the first 
                reflection and timegate the impulse response of the calibrated
                signal by. This input is optional if you want to timegate. If
                not wanting to timegate the IR, leave tgate = 0 which is the default.
                If tgate is nonzero, the IR will be gated according the input time.       
    wiener:     Boolean {True or False}; optional;
                False (default) for using direct deconvolution instead of Wiener
                deconvolution in frequency domain. If (True), the Wiener 
                deconvolution is performed. Wiener deconvolution acts as a 
                regularization which helps prevent dividing by zero allowing for 
                a more robust deconvolution while maintaining an account for any 
                system response. 
    domain:     string, Optional;
                Choice of domain performs the inverse filter in the initial step
                in either the temporal domain ('t' or 'time' or 'temporal') or 
                in the frequency domain (default) ('f' or 'freq' or 'frequency') 
                which is equivalent to determining the the cross-spectral density 
                and the auto-spectral density for the use in the deconvolution. 
                The end deconvolution always occurs by division in frequency domain.            
                
    Returns
    -------
    ht:         ndarray of float;
                Real valued impulse response (IR) of the measurement chain neglecting 
                effects of the water and tank environment through timegating only 
                direct signal with a small propagation assumption. 
    t:          ndarray of float;
                Time array for the IR h(t) in seconds (s)
    Hf:         ndarray of complex; 
                Complex two-sided Frequency Response to account for all transducer,
                amplifier, etc. in the measurement chain.     
    f:          ndarray of float
                frequency array in (Hz)
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Measurements should be taken with source and receiver close together and 
    in the center of the tank so it is easy to time gate the signal. Not too 
    close that nonlinear effects occur. The callibration position should already
    be hard coded into ESAU but can be changed manually. 
    
    The IR h(t) is determined first by timegating the signal for only direct sound
    then following Gemba(2014) eq. 3.3.6 and Farina 2000,2007 with the use of an 
    inverse filter and scaling similar to a matched filter and deconvolving in 
    order to obtain the pure delay of h(t).
    Cite: 
    "Characterization of underwater acoustic sources recorded in reverberant
    environments with application to scuba signatures" Gemba (2014) Dissertation
    
    Farina (2000)
    Farina (2007)    
    
    Often a single-sided response is desired. We find the s-sResponse
    as follows below:
        Hss = 2*Hf[0:(int(len(Hf)/2))]      #convert to single-sided FRF
        fss = f[0:int(len(f)/2)]            #convert to single-sided 
    
    last modified 5/17/2021
    """
    import numpy as np
    from ESAUResponse import IR
    
    ht = IR(cal,gen,fs,wiener=wiener,domain=domain)    #IR through deconvolution
    t = np.linspace(0,len(ht)/fs,len(ht))           #time array for ht
    
    if tgate !=0:
        print('Timegating the IR of the signal...')
        import TimeGate_UnderwaterTank as tg
        ht = tg.gatefunc(ht,fs,tgate,tb4=0.1) #cut off wall reflections
    
    #calculate the FRF from the IR and obtain the associated freq array
    print('calculating the 2-sided Frequency Response...')
    #Report the double-sided time-gated FRF of the input IR
    Hf = np.fft.fft(ht)               
    #Report the double-sided associated freq array
    f = np.fft.fftfreq(len(ht),d=1/fs)
    
    return ht,t,Hf,f









def TankResponse(rec,gen,fs,sysIR,wiener=True,domain='f'):
    """
    Parameters
    ----------
    rec:        ndarray of float;
                Received signal. Should be real valued.
    gen:        ndarray of float;
                Pure generated signal. Should be real valued.
    fs:         float; 
                Sampling frequency (Hz)
    sysIR:      ndarray;
                This is the system impulse response h(t) of the whole measurment 
                chain found between two close points using SystemResponse func.         
    wiener:     Boolean {True or False}; optional;
                False (default) for using direct deconvolution instead of Wiener
                deconvolution in frequency domain. If (True), the Wiener 
                deconvolution is performed. Wiener deconvolution acts as a 
                regularization which helps prevent dividing by zero allowing for 
                a more robust deconvolution while maintaining an account for any 
                tank response effects. 
    domain:     string, Optional;
                Choice of domain performs the inverse filter in the initial step
                in either the temporal domain ('t' or 'time' or 'temporal') or 
                in the frequency domain (default) ('f' or 'freq' or 'frequency') 
                which is equivalent to determining the the cross-spectral density 
                and the auto-spectral density for the use in the deconvolution. 
                The end deconvolution always occurs by division in frequency domain.
                 
    Returns
    -------
    H_tank:     ndarray of float; 
                Complex two-sided Greens function of Frequency Response of the Tank 
                envrionment 
    f:          ndarray of float;
                Two-sided frequency array matching the frequency response H_tank             
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    This greens function is relative to the individual positions of the Source
    and Receiver in the tank. To use this, you will need to also run the 
    SystemResponse functin to obtain the frequency response of the 
    measurement chain (transducers,etc.)
    
    last modified 5/17/2021
    """
    import numpy as np
    from ESAUResponse import IR 
    #obtain the IR of the recorded signal relative to the generated signal
    ht = IR(rec,gen,fs,wiener=wiener,domain=domain)
    #if ht.size != sysIR.size, zeropadding is necessary at the end of the smaller
    #so they are both of the same length, thus interpolating the FRF which allows
    #the components of the deconvolution to be the same size. 
    if len(ht)<len(sysIR):
        #number of zeros needed for padding to obtain same size for ht
        nzeros =int(np.abs(len(ht)-len(sysIR)))
        h = np.zeros(len(ht))
        fin = int(.999*len(ht))
        damp = int(0.0005*len(ht))
        h[0:fin] = ht[0:fin] #replace up to gate with original
        #apply half-hanning window to last portion of the array this allows for
        #the signal to more gradually ramp down to zeros to be padded. 
        h[fin:fin+damp] = ht[fin:fin+damp]*np.hanning(damp)
        #repopulate first half of that damping data keeping original array information
        h[fin:int(fin+damp/2)] = ht[fin:int(fin+damp/2)]
        #pad the end of the array with zeros making up for the difference
        ht0 = np.pad(h,(0,nzeros),'constant',constant_values=(0,0))
        sys = sysIR

    if len(sysIR)<len(ht):
        nzeros =int(np.abs(len(ht)-len(sysIR)))
        s = np.zeros(len(sysIR))
        fin = int(.999*len(sysIR))
        damp = int(0.0005*len(sysIR))
        s[0:fin] = sysIR[0:fin] #replace up to gate with original
        #apply hanning window to last portion of the array this allows for
        #the signal to more gradually ramp down to zeros to be padded. 
        s[fin:fin+damp] = sysIR[fin:fin+damp]*np.hanning(damp)
        #repopulate first half of that damping data keeping original array information
        s[fin:int(fin+damp/2)] = sysIR[fin:int(fin+damp/2)]
        sys = np.pad(s,(0,nzeros),'constant',constant_values=(0,0))
        ht0 = ht
    
    if len(sysIR) == len(ht):
        sys = sysIR
        ht0 = ht
    
    #Obtain frequency response of both the sysIR and ht for deconvolution in freq.
    Sys = np.fft.fft(sys)
    Hf = np.fft.fft(ht0)
    f = np.fft.fftfreq(len(ht0),d=1/fs) #associated frequency array        
            
    if wiener == True:
        print('performing deconvolution via Wiener deconvolution preventing dividing by zero')
        #Wiener Deconvolution deals with the near zero values which cause processing noise
        #and high frequency aliasing.
        lamb = 0.005 #scaling parameter arbitrarily chosen 
        sigma = lamb*np.mean(np.abs(Sys)) 
        WDeconv = np.conj(Sys)*Hf/(np.abs(Sys)**2+sigma**2)
        Deconv = WDeconv
    else: 
        print('Performing deconvolution via direct division in frequency domain')
        Deconv = Hf/Sys #standard deconvolution is division in freq domain
        
    Htank = Deconv
    
    return Htank,f  





