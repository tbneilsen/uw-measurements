# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:39:03 2020

@author: cvongsaw
"""

def SpectrumCalc(x,fs,pref = 20e-6, unitflag = 1): 
    import numpy as np
    """
    Compute spectrum for entire waveform x (AutoSpectral Density or Autospectrum)
    Returns chosen type of output (Autospectrum or Autospectral Density). 
    Assumes returns Autospectrum (unitflag = 1)

    Parameters
    ----------
    x:      ndarray of float;
            Signal in time domain.
    fs:     int;
            Sampling frequency. 
    pref:   float, optional;
            Reference pressure. If not given, assumed 20 microPascals as per standard. 
            Underwater standard is 1 microPascal.
    unitflag:   {0 or 1}; optional;
                0 is for autospectral density, 1 is for autospectrum.
                Default is 1 for autospectrum. 

    Returns
    -------
    Xss:    ndarray of float; 
            Single-sided complex and scaled Fourier Transform 
    Gxx:    ndarray of float;
            Single-sided real autospectrum or autospectral density
    f:      ndarray of float;
            The frequency array that matches with x
            
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    last modified 7/13/2020
    
    This code comes from Traci Neilsen as a MatLab script written 
    in order to write an autospec that doesn't average as badly as byu.autospec
    and then adapted for python. There is no windowing done within this function.
    And it does not have the time averaging.
    """  
    
    #number of samples
    ns = len(x) 
    #bin spacing(freq spacing)
    df = fs/ns  
    #frequency array
    f = np.arange(0,fs/2,df) 
    
    """
    #maximum frequency for plotting
    fmax = np.max(f)
    #determine the index of the max frequency for plotting. 
    nfmax = np.argmin(np.abs(f-fmax)) #doesnt work yet
    #time record
    t = np.linspace(0,len(x)-1)/fs
    """
 
    #Enforce zero-mean
    x = x-np.mean(x)
    
    #unscaled fourier transform to convert to frequency domain
    X = (np.fft.fft(x)) 

    #Takes first ns/2 points to make it single-sided.
    Xss = X[:int(ns/2)]
    #SPECTRUM SCALING: 0 for spectral density
    Scale = 2/ns/fs
    
    #scaled fourier transform
    Xss = Xss * np.sqrt(Scale)*np.sqrt(df)**unitflag
    #Fourier Transform turned into autospectrum
    Gxx = np.conj(Xss)* Xss 
    return Xss,Gxx,f





def CrossCalc(x,y,fs):
    import numpy as np
    """
    Compute cross spectrum and cross spectral density of two waveforms x and y
    Returns two-sided Cross Spectral Density Sab and single-sided Cross 
    Spectra Gab and the Cross Correlation

    Parameters
    ----------
    x:      ndarray of float;
            Signal in time domain. Often the received signal
    y:      ndarray of float;
            Signal in time domain. Often the generated signal
    fs:     int;
            Sampling frequency. 


    Returns
    -------
    Sab:    ndarray of float; 
            Two-sided complex cross-spectral density
    Gab:    ndarray of float;
            Single-sided real autospectrum or autospectral density
    f:      ndarray of float;
            The frequency array that matches with x
    Rab:    ndarray of float;
            Cross-Correlation
    tau:    ndarray of float;
            Time shift
            
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    last modified 7/28/2020
    
    """  
    
    #estimate from notes 560
    #Sab = a*(f)b(f)
    #a*(f) is the conjugate of the two sided FFT of a(t)
    #b(f) is the two sided FFT of b(t)
    #Cross Correlation
    #Rab(tau) = 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




"""  
  
    #AUTOSPECTRAL DENSITY
        Xss = Xss * np.sqrt(Scale)*np.sqrt(df)**outdata
        Gxx = np.mean(np.conj(Xss)*Xss,1) #Units are Pa^2/Hz
        PSD = 10*log10(np.conj(np.transpose(Gxx))/pref**2)
        return Gxx,f
         
        
        subplot(2,1,1)
        plot(t,x)
        title(titleStr)%['Pad ',num2str(padNum(jID)),' Channel ',num2str(CHnum)]);%,newline,'Waveform'])
        #xlim([min(Tw) max(Tw)])
        ylabel('Signal')
        #colorbar
        grid on
         
         
        subplot(2,1,2)
        PSD_plot = PSD(1:nfmax,:);
        h1 = semilogx(f(1:nfmax),PSD_plot);
        ylabel('Spectrum')
        xlabel('Frequency [Hz]')
        
        grid on
         
        figuresize(9,9,'inches')
        #saveFigure([titleStr, '_pt_PSD']);
"""

    
    
   







