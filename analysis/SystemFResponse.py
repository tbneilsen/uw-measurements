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

def SystemResponse(cal,gen,fs,AEgirPose,RanPose,tstart=0.5,D=0.2,T=16.0,S=0.03,pressure = True):
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
        
    last modified 02/10/2021
    """
    #TimeGate_UnderwaterTank.py timegateTank
    from TimeGate_UnderwaterTank import timegateTank
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.signal as sci 
    print('')
    print('printing calibration position details...')
    print('')
    tgate,tside,dis,c = timegateTank(AEgirPose,RanPose,D,T,S,Coordinate='tank') 
    
    #convert leading zeros time length to samples of leading zeros
    #take only majority % of the zeros allowing to view just before the signal 
    #rises with % assumed roughly by fs values for proper allowance of signal.
    #but realistically these values are mostly arbitrarily trying to approach
    #very close to the initial signal considering fs.
    if fs <= 5e5:
        percent = 0.99999
    """
    #?????????????This should be based off of the frequency band not fs????????
    """
    if fs > 5e5 and fs <= 1e6:
        percent = 0.999
    if fs > 1e6:
        percent = 0.99
    
    """
    Create a true false boolian logic to true, do the time gate, or false, dont
    """
    tstart = tstart *percent 
    start = int(tstart*fs)
    fin = int((tgate+tstart)*fs)
    
    """
    plt.figure()
    plt.plot(gen)
    plt.plot(cal-np.mean(cal))
    plt.title(f'gen and cal')
    plt.legend(['gen','cal'])
    """
    
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
    #acorr = acorr - np.mean(acorr)
    #xcorr = xcorr - np.mean(xcorr)
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
    #single sided
    Xcorr = 2*Xcorr[0:(int(len(Xcorr)/2))]
    Acorr = 2*Acorr[0:(int(len(Acorr)/2))]
    fss = f[0:int(len(f)/2)]/1000   #convert from Hz to kHz
    
    plt.figure()
    plt.plot(t,xcorr)
    plt.plot(t,acorr)
    plt.legend([r'r(t)$ \ast \overline{g(-t)}$',r'g(t)$ \ast \overline{g(-t)}$'])
    plt.xlabel('Time (ms)')
    plt.title('Components of Cal & Gen Convolution')
 
    
    #################################################################
    #********????????????????????? TO DO ***********************
    #determine Frequency band of gen from where Acorr approaches zero
    #Truncate Xcorr and Acorr for that value
    #can also do this for low frequencies to apply bandpass filter
    """
    Create a true false boolian logic to true, do the bandpass, or false, dont
    """
    cond = 0.0001*np.abs(np.max(Xcorr))
    band = np.argwhere(np.abs(Xcorr)>cond)
    hpass = min(band) #high pass frequency array index
    lpass = max(band) #low pass frequency array index
    fhigh = f[hpass] #high pass frequency
    flow = f[lpass] #low pass frequency
    
    print(f'Frequency Band = {fhigh} < f < {flow}')
    
    #deconvolution interval
    Xd = Xcorr[int(hpass):int(lpass)]
    Ad = Acorr[int(hpass):int(lpass)]
    Deconv = Xd/Ad #standard deconvolution is division in freq domain
    
    
    Acorr_dB = np.abs(10*np.log10(Acorr[int(hpass):int(lpass)])) #convert to dB for plotting visualization
    Xcorr_dB = np.abs(10*np.log10(Xcorr[int(hpass):int(lpass)])) #convert to dB for plotting visualization
    
    plt.figure()
    plt.plot(fss[int(hpass):int(lpass)],Xcorr_dB)
    plt.plot(fss[int(hpass):int(lpass)],Acorr_dB)
    plt.legend([r'fft(r(t)$ \ast \overline{g(-t)}$)',r'fft(g(t)$ \ast \overline{g(-t)}$)'])
    plt.xlabel('Frequency (kHz)')
    plt.ylabel(r'dB re 1$ \mu$Pa')
    plt.title(f'Bandpass filtered Response for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
    
    #Wiener Deconvolution (if there are noise stuggle in the deconvolution this adapts for high freq noise)
    lamb = 0.005 #scaling parameter
    #expectation or noise or SNR (Dr. Anderson uses mean instead of max and lamb = 0.9)
    #https://drive.google.com/file/d/12mprsdEDttCcWy7GOT1ySOfyf58jZCWL/view for paper by Dr. Anderson
    #using this in the deconvolution and discussing optimizing lamb
    sigma = lamb*np.max(np.abs(Ad)) 
    WDeconv = np.conj(Ad)*Xd/(np.abs(Ad)**2+sigma**2)
    
    ht = np.fft.ifft(Deconv)
    t = np.linspace(0,len(ht)/fs,len(ht))
    #convert time in sec to time in ms
    t = t*1000
    
    #################################################################
    #********????????????????????? TO DO ***********************
    #TIMEGATE ht to obtain impulse of only the direct sound before obtaining Hf
    #it is then important to zero pad the time gated signal so it is the same 
    #length as it started as again. If this is not done, then Hf has only a 
    #few frequency bins which is not enough information. 
    start = 0
    fin = int(tgate*fs*percent)
    ht1 = ht[start:fin]
    t1 = t
    #print('fin=',fin)
    #pad ht with zeros to maintain detail for Hf
    pad0 = start
    pad1 = len(ht)-len(ht1)-start
    ht1 = np.pad(ht1,(pad0,pad1),'constant',constant_values=(0))
    
    #print('ht1=',len(ht))
    #print('t1=',len(t))
    #print('ht1=',len(ht1))
    #print('t1=',len(t1))
    
    
    plt.figure()
    plt.plot(t,ht,linewidth = 4)
    plt.plot(t,ht1,linewidth = 4)
    plt.title(f'System IR for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
    plt.xlabel('Time (ms)')
    plt.legend(['non gated','gated'])
    
    Hf = np.fft.fft(ht)
    f = np.fft.fftfreq(len(ht),d=1/fs)    
    Hss = 2*Hf[0:(int(len(Hf)/2))]
    fss = f[0:int(len(f)/2)]/1000 #conver from Hz to kHz
    Hss_dB = 10*np.log10(Hss) #convert to dB for plotting visualization
    
    Hf1 = np.fft.fft(ht1)
    f1 = np.fft.fftfreq(len(ht1),d=1/fs)    
    Hss1 = 2*Hf1[0:(int(len(Hf1)/2))]
    fss1 = f1[0:int(len(f1)/2)]/1000 #conver from Hz to kHz
    Hss1_dB = 10*np.log10(Hss1) #convert to dB for plotting visualization
    
    plt.figure()
    plt.plot(fss,np.abs(Hss_dB),linewidth = 4)
    plt.plot(fss1,np.abs(Hss1_dB),linewidth = 4)
    plt.title(f'System FRF for {int(fhigh/1e3)}<f<{int(flow/1e3)} kHz Chirp')
    plt.xlabel('Frequency (kHz)')
    plt.ylabel(r'dB re 1$ \mu$Pa')
    plt.legend(['non gated','gated'])
    
    
    
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

"""
import matplotlib.pyplot as plt
f = [4.0,5.0,6.3,8.1,10.0,12.5,16.0,20.0,25.0,28.0,31.5,35.5,40.1,45.1,50.0,56.1,63.0,71.0,80.0,90.0,100.0,112.0,125.1,140.0,160.0,180.0,200.1]
Sen = [-210.0,-211.3,-211.0,-211.3,-211.3,-212.1,-212.4,-212.4,-212.0,-212.1,-212.6,-212.8,-212.3,-212.8,-213.2,-212.9,-213.4,-213.3,-213.0,-212.8,-212.1,-210.6,-211.2,-213.3,-217.3,-222.0,-220.8]
plt.figure()
plt.semilogx(f,Sen)
plt.title('BK 8103 Sensitivity')
plt.xlabel('Frequency (kHz)')
plt.ylabel(r'Sensitivity (dB re 1 V/$\mu$Pa)')
plt.xlim(4,201)
plt.minorticks_on()
plt.grid(True,which="minor",ls="-",lw=1,c='grey')
plt.grid(True,which="major",ls="-",lw=1,c='black')
plt.tick_params(which='both',width=1,color='grey')
plt.tick_params(which='major',length=7)
plt.tick_params(which='minor',length=5)
"""