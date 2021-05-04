# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:39:55 2021

@author: cvongsaw
"""

def IR(rec,gen,fs,bandpass=True,wiener=False,domain='t'):
    """
    Parameters
    ----------
    rec:        ndarray of float of size 1;
                timedomain of the received signal in pressure (if pressure = True) 
                or Voltage (if pressure = False). 
    gen:        ndarray of float;
                time domain of the generated signal in pressure (if pressure = True) 
                or Voltage (if pressure = False). 
    fs:         float; 
                Sampling frequency in Hz
    bandpass:   Boolean {True or False}; Optional;
                Default (True) for providing a bandpass filter for the deconvolution. 
                (False) for not performing a bandpass filter and taking all the noise.
    wiener:     Boolean {True or False}; optional;
                Default (False) for using direct deconvolution instead of Wiener
                deconvolution. Wiener deconvolution helps prevents the potential
                of dividing by zero allowing for a more robust deconvolution
                while maintaining any system response desired for accounting. 
                If (True), then Wiener deconvolution is performed. If bandpass
                is True, wiener is forced to be true in order to not divide by 
                zeros in the deconvolution
    domain:     string, Optional;
                Choice of domain to perform initial step of convolution where 
                an temporal inversion is applied. In the frequency domain this 
                is simply turning the deconvolution into a ratio of the cross-
                spectral density to the auto-spectral density. Defaults to 't' ('time')
                for time domain where scipy.correlate is used for the time domain.
                Choose 'f' or 'freq' or 'frequency' to perform purely in frequency 
                domain.
                
    Returns
    -------
    ht:         ndarray of float;
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
        
    last modified 5/4/2021
    """
    import numpy as np
    import scipy.signal as sci 
    
    """This choice allows for initial convolution in time domain"""
    """Deconvolution still happens only in frequency domain"""
    if domain == 'time' or 't':    
        #Gemba (2014) dissertation eq 3.3.3, 3.3.5 and 3.3.6
        #Get the temporal inverse of the generated signal so the xcorr is 
        #basically a convolution in the time domain. The temporal inversion also 
        #makes sure there is a phase inversion and thus the xcorr results in a 
        #pure delay of h(t)
        #*************************************************************************
        #rec(t)*conj(gen(-t)) = h(t)*gen(t)*conj(gen(-t)) eq 3.3.5 solve for h(t)
        gen_flip = np.conj(np.flip(gen))
        #xcorr of the recorded and temporal inverse generated (xcorr takes conj of
        #the gen_flip as part of its process) acorr is autocorrelation of the gen
        acorr = sci.correlate(gen,gen_flip,mode='full',method='auto')
        xcorr = sci.correlate(rec,gen_flip,mode='full',method='auto')
        
        #It is necessary to take these two parts into the frequency domain in order
        #to divide them out to solve for h(t). Division in the frequency domain is a
        #deconvolution in the time domain. which gives H(f) and then ifft(H(f))=h(t)
        #Convert to Frequency Domain using FFT
        """#******edit normalization of FFT???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
        Xcorr = np.fft.fft(xcorr) #double-sided frequency response
        Acorr = np.fft.fft(acorr) #double-sided frequency response
        f = np.fft.fftfreq(len(xcorr),d=1/fs)
        
    if domain == 'frequency' or 'freq' or 'f':
        """ COMPUT ALL IN FREQ DOMAIN instead of above"""
        """ Initial investigation looks like the two methods give same results"""
        Xcorr = np.fft.fft(rec)*np.conj(np.fft.fft(gen))
        Acorr = np.fft.fft(gen)*np.conj(np.fft.fft(gen))
        f = np.fft.fftfreq(len(gen),d=1/fs)
    
    
    #Single Sided 
    Xcorr = 2*Xcorr[0:(int(len(Xcorr)/2))] 
    Xd = Xcorr #renaming to call for the deconvolution later 
    Acorr = 2*Acorr[0:(int(len(Acorr)/2))]
    
    if bandpass ==True:
        #Determine Frequency band of the gen signal to apply a bandpass filter to 
        #Acorr so as to not divide by various close to zero values causing significant
        #high frequency processing noise. This cuts the leading and trailing zeros out
        #to not obtain noise outside of the signal when performing the Deconvolution.
        #**HOWEVER**, it may also remove important data from the system response. 
        """# TO DO TO DO TO DO TO DO TO DO TO DO  not perfect for sure but decent"""
        #Determine the condition in which the signal strength is small compared to max
        cond = 0.17*np.abs(np.max(Acorr))     #arbitrary coeff.
        band = np.argwhere(np.abs(Acorr)>cond)  #index values meeting above condition
        hpass = min(band)                       #high pass index of response
        lpass = max(band)                       #low pass index of response
        fhigh = f[hpass]    #high pass frequency
        flow = f[lpass]     #low pass frequency
        print(f'Frequency Band = {fhigh} < f < {flow}')
        
        #Bandpass filter only Acorr since both Gen and Genflip should not have freq. content outside bandpass filter 
        #pad Acorr with zeros to maintain original length in relation to Xcorr
        Acorr = Acorr[int(hpass):int(lpass)] #bandpass reduces size of array
        pads = hpass #start pad index up to start of freq band
        padf = len(Xcorr)-lpass #finish pad index after end of band
        Acorr = np.pad(Acorr,(int(pads),int(padf)),'constant',constant_values=(0)) #regain original size with zero padding
        Ad = Acorr #renaming to call for the deconvolution later 
        
        #Since the bandpass puts a lot of zeros to be divided in the deconvolution, 
        #we will set wiener = True such that we do not divide by zero****
        wiener = True
        
    else: 
        #When not bandpassing the signal (which takes out potential effects from 
        #the system response) repass Xd and Ad as the original arrays
        Xd = Xcorr 
        Ad = Acorr
    
    if wiener == True:
        print('Performing deconvolution via Wiener deconvolution preventing dividing by zero')
        #Wiener Deconvolution deals with the near zero values which cause processing noise.
        #The noise is often high frequency noise. 
        #This can be considered an alternate option to avoiding dividing by zero other than
        #using the bandpass method which loses potential system noise you may which to 
        #account for.
        lamb = 0.005 #scaling parameter arbitrarily chosen for now
        #expectation or noise or SNR (Dr. Anderson uses mean instead of max and lamb = 0.9)
        #https://drive.google.com/file/d/12mprsdEDttCcWy7GOT1ySOfyf58jZCWL/view for paper 
        #by Dr. Anderson using this in the deconvolution and discussing optimizing lamb. 
        #Wiener deconvolution for image restoration also typically uses the mean.
        sigma = lamb*np.mean(np.abs(Ad)) #expectation or noise or SNR
        WDeconv = np.conj(Ad)*Xd/(np.abs(Ad)**2+sigma**2)
        Deconv = WDeconv
    else: 
        print('Performing deconvolution via direct division in frequency domain')
        #Perform standard deconvolution by direct division in frequency domain. 
        Deconv = Xd/Ad
        
    #bring back to time domain with inverse fast fourier transform (IFFT)
    """#******edit normalization of FFT???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    ht = np.fft.ifft(Deconv) 
    
    return ht
    
    





    
def SysResponse(cal,gen,fs,AEgirPose,RanPose,tgate=0,bandpass=True,wiener=False):
    """
    Parameters
    ----------
    cal:        ndarray of float;
                Received calibration signal 
    gen:        ndarray of float;
                Pure generated signal
    fs:         float; 
                Sampling frequency (Hz)
    AEgir_Pose: tuple; 
                AEgir TCP position (x,y,z) in the 'tank' frame 
    Ran_pose:   tuple;
                Ran TCP position (x,y,z) in the 'tank' frame
    tgate:      float,Optional;
                Time of the first wall reflection determined through timegateTank. 
                This is the time we will use to determine the impule of the 1st 
                reflection and timegate the impulse response of the calibrated
                signal by. This input is optional if you want to timegate. If
                not wanting to timegate the IR, leave tgate = 0 which is the default.
                If tgate is nonzero, the IR will be gated according the input time.       
    bandpass:   Boolean {True or False}; optional;
                Default (True) for providing a bandpass filter for the deconvolution. 
                (False) for not performing a bandpass filter and taking all the information.
    wiener:     Boolean {True or False}; optional;
                Default (False) for using direct deconvolution instead of Wiener
                deconvolution. Wiener deconvolution helps prevents the potential
                of dividing by zero allowing for a more robust deconvolution
                while maintaining any system response desired for accounting. 
                If (True), then Wiener deconvolution is performed. 
                
    Returns
    -------
    ht:         ndarray of float;
                Real impulse response (IR) of the measurement chain neglecting 
                effects of the water and tank environment through timegating only 
                direct signal with small propagation assumption. 
    t:          ndarray of float;
                Time array for the IR h(t) in (ms)
    Hf:         ndarray of complex; 
                Complex two-sided Frequency Response to account for all transducer,
                amplifier, etc. in the measurement chain.     
    f:          ndarray of float
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
        
    last modified 5/4/2021
    """
    import numpy as np
    from ESAUResponse import IR
    
    ht = IR(cal,gen,fs,bandpass=bandpass,wiener=wiener) #IR through deconvolution
    t = np.linspace(0,len(ht)/fs,len(ht))               #time array for ht
    t = t*1000                                          #convert time in sec to time in ms
    
    if tgate !=0:
        print('Timegating the IR of the signal...')
        #take only majority % of the zeros allowing to view just before the signal 
        #rises with % assumed roughly by fs values for an allowance of the signal.
        #but realistically these values are mostly arbitrarily trying to approach
        #very close to the initial signal considering fs.
        """#?????????This should be based off of the frequency band not fs!!!!"""
        
        """
        Might be good to look into constant fraction descrimination to find peak
        around the tgate time"""
        
        import TimeGate_UnderwaterTank as tg
        ht1 = tg.gatefunc(ht,fs,tgate,tb4=0.1)
        
        """
        tb4gate = 0.1 #ms before first reflection tgate
        Nb4gate = tb4gate/1000 *fs #convert to samples before gating
        
        fin = int(tgate*fs-Nb4gate)
        start = 0 #start the time gating allowing everything from the beginning of ht
        #convert time length to samples to determine the finish cutoff of ht
        #fin = int(tgate*fs*percent)
        #cut off the IR before the first reflection being index "fin"
        ht1 = np.zeros(len(ht))
        ht1[start:fin] = ht[start:fin] #replace up to gate with original
        damp = int(0.05/1000*fs) #0.05ms of damping converted to samples
        #apply hanning window to portion of the array following original cutoff
        #this allows for the signal to more gradually ramp down to zeros. 
        ht1[fin:fin+damp] = ht[fin:fin+damp]*np.hanning(damp)
        #repopulate first half of that damping data keeping original array information
        ht1[fin:int(fin+damp/2)] = ht[fin:int(fin+damp/2)]
        """
        
    
    #calculate the FRF from the IR and obtain the associated freq array
    print('calculating the 2-sided Frequency Response...')
    #report the double-sided time-gated FRF
    Hf = np.fft.fft(ht1)               #FRF of the time-gated system IR
    #report the double-sided associated freq array
    f = np.fft.fftfreq(len(ht),d=1/fs)#associated frequency array     
    
    #Often a single-sided response is desired. We find the s-sResponse (unused)
    #as follows below:
        #Hss = 2*Hf[0:(int(len(Hf)/2))]   #convert to single-sided FRF
        #fss = f[0:int(len(f)/2)]/1000    #convert to single-sided from Hz to kHz 
    
    #reset variables to be ready for output. 
    ht = ht1    #report the time-gated IR
    
    return ht,t,Hf,f









def TankResponse(rec,gen,fs,sysIR,bandpass=True,wiener=True):
    #!!!!!bandpass should default to False and wiener to True for TankResponse???
    """
    Parameters
    ----------
    rec:        ndarray of float;
                Received signal  
    gen:        ndarray of float;
                Pure generated signal
    fs:         float; 
                Sampling frequency
    sysIR:      ndarray;
                This is the system impulse response of the whole measurment 
                chain found between two close points using 
                SystemFResponse.py SystemResponse         
    bandpass:   Boolean {True or False}; optional;
                Default (True) for providing a bandpass filter for the deconvolution. 
                (False) for not performing a bandpass filter and taking all the information.
    wiener:     Boolean {True or False}; optional;
                Default (False) for using direct deconvolution instead of Wiener
                deconvolution. Wiener deconvolution helps prevents the potential
                of dividing by zero allowing for a more robust deconvolution
                while maintaining any system response desired for accounting. 
                If (True), then Wiener deconvolution is performed. If bandpass 
                is true, then wiener deconvolution is forced to happen for the 
                calculation of the IR of the recorded signal but not for the 
                deconvolution of the system response. 
                 
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
    
    This greens function is obviously relative to the individual points of
    Source and Receiver in the tank. To use this, you will need to also run
    the function SystemResponse to obtain the frequency response of the 
    measurement chain (transducers,etc.)
    
    last modified 4/12/2021
    """
    import numpy as np
    from ESAUResponse import IR 
    
    #obtain the IR of the recorded signal relative to the generated signal
    ht = IR(rec,gen,fs,bandpass=bandpass,wiener=wiener)
    #if one response is smaller than the other, we need to zeropad the end so 
    #they are both of the same length, thus interpolating the FRF which allows
    #the components of the deconvolution to be the same size. 
    if len(ht)<len(sysIR):
        #number of zeros needed for padding to obtain same size
        nzeros =int(np.abs(len(ht)-len(sysIR)))
        
        #apply hanning window to last portion of the array this allows for
        #the signal to more gradually ramp down to zeros to be padded. 
        h = np.zeros(len(ht))
        fin = int(.999*len(ht))
        damp = int(0.0005*len(ht))
        h[0:fin] = ht[0:fin] #replace up to gate with original
        h[fin:fin+damp] = ht[fin:fin+damp]*np.hanning(damp)
        #repopulate first half of that damping data keeping original array information
        h[fin:int(fin+damp/2)] = ht[fin:int(fin+damp/2)]
        
        ht = h
        ht = np.pad(ht,(0,nzeros),'constant',constant_values=(0,0))
        sys = sysIR
        Sys = np.fft.fft(sys)
    #Obtain the frequency response of the recorded measurement by converting IR 
    #to Frequency domain in order to later be deconvolved by the system Response
    Hf = np.fft.fft(ht)
    f = np.fft.fftfreq(len(ht),d=1/fs) #associated frequency array
    
    #Need to interpolate SystemResponse to match size of rec & gen lined up with 
    #the given frequency array. This is especially necessary where cal measure and
    #recorded measure are different signals. To do this, convert back to time 
    #domain to avoid interpolation of complex values using ifft. Pad the time 
    #domain with zeros on the end in order get back to the correct size and return 
    #the freq. domain using fft. 
    if len(sysIR) < len(ht):
        htsys = sysIR 
        nzeros =int(np.abs(len(ht)-len(htsys)))
        
        #apply hanning window to last portion of the array this allows for
        #the signal to more gradually ramp down to zeros to be padded. 
        s = np.zeros(len(sys))
        fin = int(.999*len(sys))
        damp = int(0.0005*len(sys))
        s[0:fin] = sys[0:fin] #replace up to gate with original
        s[fin:fin+damp] = sys[fin:fin+damp]*np.hanning(damp)
        #repopulate first half of that damping data keeping original array information
        s[fin:int(fin+damp/2)] = sys[fin:int(fin+damp/2)]
        sys = s
        sys = np.pad(htsys,(0,nzeros),'constant',constant_values=(0,0))
    
        Sys = np.fft.fft(sys)
    
    if len(sysIR) == len(ht):
        Sys = np.fft.fft(sysIR)
    #Greens function or frequency response of the tank at specific AEgir & Ran
    #positions. eq. 1.7.1 Leishman 560 notes 2019.
    #Htank = Rec/(Gen*Sys) #original idea directly from notes dividing out the 
    #system response from the recorded response to have a more pure response
    #aka a deconvolution in the time domain. Adapt for potential of dividing by 
    #zeros when using complex numbers by multiplying by 1 as (conj(sys)/conj(sys))
    #and (conj(Gen)/conj(Gen)) obtaining the magnitude on bottom. 

    if wiener == True:
        print('performing deconvolution via Wiener deconvolution preventing dividing by zero')
        #Wiener Deconvolution (if there are noise stuggle in the deconvolution this adapts for 
        #high freq noise). However this is limited to only linear time-invariant signals
        lamb = 0.6 #scaling parameter arbitrarily chosen
        #expectation or noise or SNR (Dr. Anderson uses mean instead of max and lamb = 0.9)
        #https://drive.google.com/file/d/12mprsdEDttCcWy7GOT1ySOfyf58jZCWL/view for paper by Dr. Anderson
        #using this in the deconvolution and discussing optimizing lamb
        #Wiener deconvolution for image restoration also typically uses the mean instead of max.
        sigma = lamb*np.mean(np.abs(Sys)) 
        WDeconv = np.conj(Sys)*Hf/(np.abs(Sys)**2+sigma**2)
        Deconv = WDeconv
    else: 
        print('Performing deconvolution via direct division in frequency domain')
        Deconv = Hf/Sys #standard deconvolution is division in freq domain
    
    Htank = Deconv
    
    return Htank,f  





