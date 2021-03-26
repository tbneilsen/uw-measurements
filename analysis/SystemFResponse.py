# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:33:19 2020

@author: cvongsaw
"""

def SystemResponse(cal,gen,fs,AEgirPose,RanPose,tstart=0.5,D=0.2,T=16.0,S=0.03,bandpass=True,gate=True,wiener=False):
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
    bandpass:   Boolean {True or False}; optional;
                Default (True) for providing a bandpass filter for part of the 
                deconvolution(Auto-Convolution of gen(t)).(False) for not 
                performing a bandpass filter and taking all the information.
    gate:       Boolean {True or False}; optional;
                Default (True) for time-gating the deconvolved IR h(t) based on 
                ray theory geometries via TimeGate_UnderwaterTank. The timegated and 
                non-timegated values are compared for IR and FRF if True. 
                If (False) then no timegating occurs.
    wiener:     Boolean {True or False}; optional;
                Default (False) for using direct deconvolution instead of Wiener
                deconvolution. Wiener deconvolution helps prevents the potential
                of dividing by zero allowing for a more robust deconvolution
                while maintaining any system response desired for accounting. 
                If (True), then Wiener deconvolution is performed. If bandpass
                is True, wiener is forced to be true in order to not divide by 
                zeros in the deconvolution
        
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
    
    last modified 03/26/2021
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.signal as sci 
    from TimeGate_UnderwaterTank import timegateTank
    import pdb
    plot = False
    
    #Gemba (2014) dissertation eq 3.3.3, 3.3.5 and 3.3.6
    #Get the temporal inverse of the generated signal so the xcorr is 
    #basically a convolution in the time domain. The temporal inversion also 
    #makes sure there is a phase inversion and thus the xcorr results in a 
    #pure delay of h(t)
    #*************************************************************************
    #rec(t)*conj(gen(-t)) = h(t)*gen(t)*conj(gen(-t)) eq 3.3.5 solve for h(t)
    gen_flip = np.flip(gen)
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
    #filtering)
    Xcorr = np.fft.fft(xcorr)
    Acorr = np.fft.fft(acorr)
    f = np.fft.fftfreq(len(xcorr),d=1/fs)
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
    """This is not perfect and does not do the best just yet at at determining band"""
    """This could just be another input of the chirp bandwidth to be more precise"""
    #determine Frequency band of gen from where Acorr approaches zero
    #Truncate Acorr for that bandwidth and pad with zeros to retain size.
    #The bandpass filter is not applied to Xcorr because it contains system 
    #noise that we wish to account for (though some of the high freq noise is
    #definitely an artifact of the convolution).  
    if bandpass == True:
        """ attempted to determine bandpass from Gen Signal to be more precise, but its the wrong size
        Gen = np.fft.fft(gen)
        Gss = Gen[0:int(len(Gen)/2)]/1000 
        cond = 0.01*(np.max(Gss))
        band = np.argwhere((Gss)>cond)
        hpass = min(band)   #high pass frequency array index
        lpass = max(band)   #low pass frequency array index
        fhigh = f[hpass]    #high pass frequency
        flow = f[lpass]     #low pass frequency
        print(f'Frequency Band = {fhigh} < f < {flow}')
        """
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
        
        #for plotting each component to check if it is functioning well. 
        Acorr_dB = np.abs(10*np.log10(Acorr[int(hpass):int(lpass)])) #bandpass & convert to dB for plotting visualization
        Acorr_dB = np.pad(Acorr_dB,(int(pads),int(padf)),'constant',constant_values=(0)) #pad to retain length
        Xcorr_dB = np.abs(10*np.log10(Xcorr[:])) #convert to dB for plotting visualization
        
        if plot == True:
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
        Acorr_dB = np.abs(10*np.log10(Acorr)) #convert to dB for plotting visualization
        Xcorr_dB = np.abs(10*np.log10(Xcorr)) #convert to dB for plotting visualization
        
        if plot == True:
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
        lamb = 0.005 #scaling parameter arbitrarily chosen
        #expectation or noise or SNR (Dr. Anderson uses mean instead of max and lamb = 0.9)
        #https://drive.google.com/file/d/12mprsdEDttCcWy7GOT1ySOfyf58jZCWL/view for paper by Dr. Anderson
        #using this in the deconvolution and discussing optimizing lamb
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
    
    ###########################################################################
    #TIMEGATE ht to obtain impulse of only the direct sound before obtaining Hf
    #it is then important to zero pad the time gated signal so it is the same 
    #length as it started as again. If this is not done, then Hf has only a 
    #few frequency bins which is not enough information. 
    print('')
    print('printing calibration position details...')
    print('')
    #Obtain earliest reflection times according to ray theory geometry of the tank
    tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank') 
    
    if gate == True:
        #take only majority % of the zeros allowing to view just before the signal 
        #rises with % assumed roughly by fs values for an allowance of the signal.
        #but realistically these values are mostly arbitrarily trying to approach
        #very close to the initial signal considering fs.
        #This also plots the timegated and non timegated IR & FRF for comparison
        """#?????????This should be based off of the frequency band not fs!!!!"""
        """1T of fmax(most dominant freq) buffer before 1st reflection"""
        #min period is determined from max frequency or frequency of interest
        Tmin = 1/fhigh 
        """temporary solution since the below adjustment is not working well and
        definitely not determined generically for any measurement. 
        Might be good to look into constant fraction descrimination to find peak
        around the tgate time"""
        if fs <=1e6:
            Nmin = 10 #fs = 1M
        else:
            Nmin = 1e5 #fs = 10M
        #Nmin = Tmin*fs #number of samples to shift by based on min period
        #shift gate before by Nmin samples
        fin = int(tgate*fs-Nmin)
        start = 0 #start the time gating allowing everything from the beginning of ht
        #convert time length to samples to determine the finish cutoff of ht
        #fin = int(tgate*fs*percent)
        #cut off the IR before the first reflection being index "fin"
        ht1 = ht[start:fin]
        #pad ht with zeros to maintain original length (maintains detail for Hf)
        pad0 = start
        pad1 = len(ht)-len(ht1)-start
        ht1 = np.pad(ht1,(pad0,pad1),'constant',constant_values=(0))
        
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
    
    
    """
    THIS IS THE OLDER CODE THAT DID NOT USE AN INVERSE FILTER TO ADJUST 
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
    """
    #direct method of deconvolution division in freq domain without Weiner adjustment
    Hf = np.conj(Gen)*Cal_mod/(np.abs(Gen)**2)
    ht = np.fft.ifft(Hf)
    """
    
    return ht,t,Hf,f