# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:33:19 2020

@author: cvongsaw
"""

def SystemResponse(cal,gen,fs,AEgirPose,RanPose,tgate=0,bandpass=True,wiener=False,domain='t'):
    """
    Parameters
    ----------
    cal:    ndarray of float;
            Received calibration signal. 
    gen:    ndarray of float;
            Pure generated signal. 
    fs:     float; 
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
                Default (True) for providing a bandpass filter for part of the 
                deconvolution(Auto-Convolution of gen(t)).(False) for not 
                performing a bandpass filter and taking all the information.
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
    ht:     ndarray of float;
            Real impulse response (IR) of the measurement chain neglecting 
            effects of the water and tank environment through timegating only 
            direct signal with small propagation assumption. 
    t:      ndarray of float;
            Time array for the IR h(t) reported in (ms)
    Hf:     ndarray of complex; 
            Complex two-sided Frequency Response to account for all transducer,
            amplifier, etc. in the measurement chain.     
    f:      ndarray of float
            Two-sided frequency array in (Hz)
    
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
    
    A bandpass filter is applied to the Auto-Convolution of gen(t) in order to
    remove high frequency content brought about as an artifact of the convolution
    
    last modified 03/29/2021
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.signal as sci 
    import pdb
    plot = False
    
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
        xcorr = sci.correlate(cal,gen_flip,mode='full',method='auto')
        #acorr = acorr - np.mean(acorr) #option to zero-mean removing DC offset
        #xcorr = xcorr - np.mean(xcorr) #option to zero-mean removing DC offset
        t = np.linspace(0,len(acorr)/fs,len(acorr))
        #convert time in sec to time in ms
        t = t*1000
        
        #convert to Frequency domain to perform the devoncolution
        #ht is the xcorr of rec w/ -t of gen divided by acorr which is a 
        #deconvolution in the time domain which solves for h(t) (this uses inverse
        #time filtering)
        Xcorr = np.fft.fft(xcorr)
        Acorr = np.fft.fft(acorr)
        f = np.fft.fftfreq(len(xcorr),d=1/fs)
    
    if domain == 'frequency' or 'freq' or 'f':
        """ COMPUT ALL IN FREQ DOMAIN instead of above"""
        """ Initial investigation looks like the two methods give same results"""
        Xcorr = np.fft.fft(cal)*np.conj(np.fft.fft(gen))
        Acorr = np.fft.fft(gen)*np.conj(np.fft.fft(gen))
        f = np.fft.fftfreq(len(gen),d=1/fs)
    
    #obtain single sided values
    Xcorr = 2*Xcorr[0:(int(len(Xcorr)/2))]
    Xd = Xcorr #renaming to call for the deconvolution later 
    Acorr = 2*Acorr[0:(int(len(Acorr)/2))]
    fss = f[0:int(len(f)/2)]/1000     #convert from Hz to kHz 
    
    if plot == True:
        plt.figure()
        plt.plot(t,xcorr)
        plt.plot(t,acorr)
        plt.legend([r'r(t)$ \ast \overline{g(-t)}$',r'g(t)$ \ast \overline{g(-t)}$'])
        plt.xlabel('Time (ms)')
        plt.title('Components of Cal & Gen Convolution')
 
    #################################################################
    #********????????????????????? TO DO ***********************
    """#This is not perfect and does not do the best just yet at at determining band"""
    """This could just be another input of the chirp bandwidth to be more precise"""
    #determine Frequency band of gen from where Acorr approaches zero
    #Truncate Acorr for that bandwidth and pad with zeros to retain size.
    #The bandpass filter is not applied to Xcorr because it contains system 
    #noise that we wish to account for (though some of the high freq noise is
    #definitely an artifact of the convolution).  
    if bandpass == True:
        
        #set condition where energy is high for frequencies of interest
        #coefficient is arbitrarily chosen to obtain close enough band
        #this can be improved by having the band input and found
        cond = 0.17*(np.max(Acorr)) #set condition where energy is high for frequencies of interest
        band = np.argwhere((Acorr)>cond) #determine idx values where the condition is true
        hpass = min(band)   #high pass frequency array index
        lpass = max(band)   #low pass frequency array index
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
        
        if plot == True:
            #for plotting each component to check if it is functioning well. 
            Acorr_dB = np.abs(10*np.log10(Acorr[int(hpass):int(lpass)])) #bandpass & convert to dB for plotting visualization
            Acorr_dB = np.pad(Acorr_dB,(int(pads),int(padf)),'constant',constant_values=(0)) #pad to retain length
            Xcorr_dB = np.abs(10*np.log10(Xcorr[:])) #convert to dB for plotting visualization
            plt.figure()
            plt.plot(fss,Xcorr_dB)
            plt.plot(fss,Acorr_dB)
            plt.legend([r'fft(r(t)$ \ast \overline{g(-t)}$)',r'fft(g(t)$ \ast \overline{g(-t)}$)'])
            plt.xlabel('Frequency (kHz)')
            plt.ylabel(r'dB re 1$ \mu$Pa')
            plt.title(f'Bandpass filtered Response for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
            
    else:
        Xd = Xcorr #renaming to call for the deconvolution later 
        Ad = Acorr #renaming to call for the deconvolution later 
        
        if plot == True:
            Acorr_dB = np.abs(10*np.log10(Acorr)) #convert to dB for plotting visualization
            Xcorr_dB = np.abs(10*np.log10(Xcorr)) #convert to dB for plotting visualization
            plt.figure()
            plt.plot(fss,Xcorr_dB)
            plt.plot(fss,Acorr_dB)
            plt.legend([r'fft(r(t)$ \ast \overline{g(-t)}$)',r'fft(g(t)$ \ast \overline{g(-t)}$)'])
            plt.xlabel('Frequency (kHz)')
            plt.ylabel(r'dB re 1$ \mu$Pa')
            plt.title(f'Response')
        
    if wiener == True:
        print('performing deconvolution via Wiener deconvolution preventing dividing by zero')
        #Wiener Deconvolution (if there are noise stuggle in the deconvolution this adapts for 
        #high freq noise). However this is limited to only linear time-invariant signals.
        lamb = 0.005 #scaling parameter arbitrarily chosen (HELPS DETERMINE how much deconvolution blows up)
        #expectation or noise or SNR (Dr. Anderson uses mean instead of max and lamb = 0.9)
        #https://drive.google.com/file/d/12mprsdEDttCcWy7GOT1ySOfyf58jZCWL/view for paper by Dr. Anderson
        #using this in the deconvolution and discussing optimizing lamb. Wiener deconvolution 
        #for image restoration also typically uses the mean.
        #sigma can also be known as the SNR
        sigma = lamb*np.mean(np.abs(Ad)) 
        WDeconv = np.conj(Ad)*Xd/(np.abs(Ad)**2+sigma**2)
        Deconv = WDeconv
    else: 
        print('Performing deconvolution via direct division in frequency domain')
        Deconv = Xd/Ad #standard deconvolution is division in freq domain
    
    ht = np.fft.ifft(Deconv)                #IR from deconvolution
    t = np.linspace(0,len(ht)/fs,len(ht))   #time array for ht
    t = t*1000                              #convert time in sec to time in ms
    
    #calculate the non-timegated FRF of the Cal measurement.
    #this is done for comparison and analysis. 
    Hf = np.fft.fft(ht)                 #convert from time domain to frequency via FFT (double-sided)
    f = np.fft.fftfreq(len(ht),d=1/fs)  #obtain the double-sided associated frequency array
    Hss = 2*Hf[0:(int(len(Hf)/2))]      #obtained the single sided FRF from the double sided
    fss = f[0:int(len(f)/2)]/1000       #convert from Hz to kHz
    Hss_dB = 10*np.log10(Hss)           #convert to dB for plotting visualization
    
    plt.figure()
    plt.plot(fss,Hss_dB)
    plt.title('Deconvolved TTS Response')
    plt.xlabel('Frequency (kHz)')
    plt.ylabel(r'Level (dB re 1 $\mu Pa$)')
    plt.grid()
    
    ###########################################################################
    #TIMEGATE ht to obtain impulse of only the direct sound before obtaining Hf
    #it is then important to zero pad the time gated signal so it is the same 
    #length as it started as again. If this is not done, then Hf has only a 
    #few frequency bins which is not enough information. 

    #Obtain earliest reflection times according to ray theory geometry of the tank
    #prior to using this code by using timegateTank() from TimeGate_UnderwaterTank
    #as follows:
    #tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank') 
    
    if tgate != 0:
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
        
        
        if plot == True:
            #plot timegated & non-timegated IR for reference and comparison. 
            plt.figure()
            plt.plot(t,ht,linewidth = 4)    #plot the non-timegated IR
            plt.plot(t,ht1,linewidth = 4)   #plot the timegated IR
            if bandpass == True:
                plt.title(f'System IR for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
            else:
                plt.title(f'System IR')
            plt.xlabel('Time (ms)')
            plt.legend(['non gated','gated'])
    
    
        #calculate the FRF from the time-gated IR and obtain the associated freq array
        Hf1 = np.fft.fft(ht1)               #FRF of the time-gated system IR
        f1 = np.fft.fftfreq(len(ht1),d=1/fs)#associated frequency array    
        Hss1 = 2*Hf1[0:(int(len(Hf1)/2))]   #convert to single-sided FRF
        fss1 = f1[0:int(len(f1)/2)]/1000    #convert to single-sided from Hz to kHz
        Hss1_dB = 10*np.log10(Hss1)         #convert to dB for plotting visualization
        
        if plot == True:
            plt.figure()
            plt.plot(fss,np.abs(Hss_dB),linewidth = 4)
            plt.plot(fss1,np.abs(Hss1_dB),linewidth = 4)
            if bandpass == True:
                plt.title(f'System FRF for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
            else:
                plt.title(f'System FRF')
            plt.xlabel('Frequency (kHz)')
            plt.ylabel(r'dB re 1$ \mu$Pa')
            plt.legend(['non gated','gated'])
        
        #reset variables to be ready for output. 
        ht = ht1    #report the time-gated IR
        Hf = Hf1    #report the double-sided time-gated FRF
        f = f1      #report the double-sided associated freq array
    
    
    return ht,t,Hf,f