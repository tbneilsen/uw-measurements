# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:25:42 2019

@author: cvong
"""
###############################################################################
##################          UW Measurements       #############################
#########     Obtain time waveforms, Auto & Cross Spec,     ###################
#########     FRF, impulse response, spectrogram            ###################
###############################################################################
##################          TC4034 receiver       #############################
##################           TC4038 Source        #############################
###############################################################################
###############################################################################
import byuarglib as byu
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import ifft

import sys
sys.path.insert(1, "C:\Users\cvong\Box\UW Research\Code\python-plotting-lib")
from fig_defaults import set_defaults as set_fig_defaults

_ = set_fig_defaults("paper")
######################################
########### Formatting ###############
#loaded_signal = byu.binfileload(path, IDname, IDnum, CHnum)
#number of files if beginning at ID001
#loads as pressures in (Pa)
######################################
files_begin     = 12     #first file (python indexes from 0)
files_end       = 13    #last file (python indexes up to last value)
path = "C:/Users/cvong/Box/UW Research/Underwater Measurements/2020-02-06 Tank test pulses 2"
fs = 1240000            #sample frequency (Hz)
t_rec_len = 10          #time record length (s)
Ts = 1/fs               #sampling period (s)
N = fs*t_rec_len        #total number of samples 
n=15                    #number of bins
ns = 2**15              #samples per block
                                        #######################################
                                        #set time record window to view...
sampleview = N                          #how many samples you want to focus...  
t_firstwindow = Ts*sampleview           #on."N" being the full value. Adjust...
time = np.arange(0,t_firstwindow,Ts)    #for viewing first window of samples. 
                                        #######################################
                                        
#load all binfiles for generated (CH0) and received (CH1) signals
all_sig_gen = np.empty((files_end-files_begin,N))
all_sig_rec = np.empty((files_end-files_begin,N))

for ii in range(files_begin,files_end):
    all_sig_gen[ii-files_begin,:] = byu.binfileload(path,"ID",ii,0)   #signals generated
    all_sig_rec[ii-files_begin,:] = byu.binfileload(path,"ID",ii,1)   #signals recorded

for ii in range(files_begin,files_end):
    plt.figure()
    plt.plot(time, all_sig_gen[ii-files_begin])
    plt.title("T-Waveform for sig. generated {}".format(ii))
    plt.figure()
    plt.plot(time, all_sig_rec[ii-files_begin])
    plt.title("T-Waveform for sig. recorded {}".format(ii))
plt.show()


#########################
#Parameters byu.autospec#
#########################
#    ----------
#    x : ndarray of float
#        Time series signal. Should be a single dimensional array.
#    fs : int
#        Sampling frequency used to collect the signal
#    ns : int, optional
#        Number of samples per block. Defaults to 2^15 if not specified. If
#        input ns is larger than the number of samples, the nearest lower
#        multiple of 2 is selected.
#    N : int, optional
#        Total number of samples to test. If N is not an integer multiple of ns,
#        the samples less than ns in the last block are discarded. Defaults to
#        lowest power of 2 if not specified.
#    unitflag : {0 or 1}, optional
#        1 for autospectrum and 0 for autospectral density. Defaults to
#        autospectral density
#########################
#    Returns
#########################
#    Gxx : ndarray of float
#        the single-sided finished autospectrum or autospectral density
#    f : ndarray of float
#        the frequency array that matches with x
#    OASPL : float
#        sound pressure level
#########################

Gxx_gen = np.empty((files_begin,files_end))
OASPL_gen = np.empty((files_begin,files_end))
Gxx_rec = np.empty((files_begin,files_end))
OASPL_rec = np.empty((files_begin,files_end))
for ii in range(files_begin,files_end):
    Gxx_gen[ii,:], f_new, OASPL_gen[ii] = byu.autospec(all_sig_gen[ii], fs, ns, N)
    Gxx_rec[ii,:], f_new, OASPL_rec[ii] = byu.autospec(all_sig_rec[ii], fs, ns, N)

#########################
#Parameters byu.crossspec
#########################
#    x : ndarray of float
#        Time series signal. Should be a single dimensional array.
#    y : ndarray of float
#        Time series signal. Should be a single dimensional array.
#    fs : int
#        Sampling frequency used to collect the signal
#    ns : int, optional
#        Number of samples per block. Defaults to 2^15 if not specified. If
#        input ns is larger than the number of samples, the nearest lower
#        multiple of 2 is selected.
#    N : int, optional
#        Total number of samples to test. If N is not an integer multiple of ns,
#        the samples less than ns in the last block are discarded. Defaults to
#        lowest power of 2 if not specified.
#    unitflag : {0 or 1}, optional
#        1 for autospectrum and 0 for autospectral density. Defaults to
#        autospectral density
#########################
#    Returns
#########################
#    Gxy : ndarray of float
#        Gxy is the single-sided finished cross spectrum or cross
#        spectral density
#    f : ndarray of float
#        f is the frequency array that matches with x
#########################

Gxy = np.empty(files_begin,files_end)
for ii in range(files_begin,files_end):
    Gxy[ii,:], f_new = byu.crossspec(all_sig_gen[ii], all_sig_rec[ii], fs, ns, N)

#########################
#Freq Response Function #
#########################
H = np.empty(files_begin,files_end)
for ii in range(files_begin,files_end):
    H[ii,:] = Gxy[ii]/Gxx_gen[ii]
#########################
# inverse fasr fourier  #
# Parameters ifft       #
#########################
#   x: array_like
#       Transformed data to invert.
#   n: int, optional
#       Length of the inverse Fourier transform. If n < x.shape[axis], 
#       x is truncated. If n > x.shape[axis], x is zero-padded. 
#       The default results in n = x.shape[axis].
#   axis: int, optional
#       Axis along which the ifft’s are computed; 
#       the default is over the last axis (i.e., axis=-1).
#   overwrite_x: bool, optional
#       True, the contents of x can be destroyed; the default is False.
#########################
#Returns
#########################
#   ifft: ndarray of floats
#       The inverse discrete Fourier transform.
#########################

ImpRes = np.empty(files_begin,files_end)
for ii in range(files_begin,files_end):
    ImpRes[ii] = ifft(H[ii])
    



#########################
#Parameters byu.specgram#
#########################
#x : ndarray of float
#    Tme series signal. Should be a single dimensional array.
#fs : int
#    Sampling frequency of the signal
#ns : int
#    Number of samples per block. Must be divisible by 4. (use 2**n)
#pct : int or float, optional
#    Percentage overlap. Defaults to 0%
#unitflag : {1 or 0}, optional
#    1 for autospectrum, 0 for autospectral density. Defaults to
#    autospectral density.
#########################
##Returns byu.specgram###
#########################
#Gxx : ndarray of float
#    Single-sided autospectrum or autospectral density for each time block,
#    depending on unitflag
#t : ndarray of float
#    t is the time array corresponding to the first sample in each block
#f : ndarray of float
#    f is the frequency array
#runOASPL : ndarray of float
#    running OASPL for each time blockarray
##########################

Gxx_gen = np.empty(files_begin,files_end)
runOASPL_gen = np.empty(files_begin,files_end)
Gxx_rec = np.empty(files_begin,files_end)
runOASPL_rec = np.empty(files_begin,files_end)

for ii in range(files_begin,files_end):
    Gxx_gen[ii],t_gen,f_gen,runOASPL_gen[ii] = byu.specgram(all_sig_gen,fs,ns)
    Gxx_rec[ii],t_rec,f_rec,runOASPL_rec[ii] = byu.specgram(all_sig_rec,fs,ns)
"""
for ii in range(files_begin,files_end):
    plt.figure()
    plt.pcolormesh(t_gen[ii-file_begins], f_gen[ii-file_begins],Gxx_gen[ii-file_begins])
    plt.title("Spectrogram for sig. generated {}".format(ii))
    plt.figure()
    plt.pcolormesh(t_rec[ii-file_begins], f_rec[ii-file_begins],Gxx_rec[ii-file_begins])
    plt.title("Spectrogram for sig. recorded {}".format(ii))
plt.show()
"""










"""
#autospec punches out the autospectrum in pressures, the frequencies and the overall SPL
Gxx1, f1, OASPL1 = byu.autospec(ID1,fs,102400)
Gxx2, f2, OASPL2 = byu.autospec(ID2,fs,102400)
Gxx3, f3, OASPL3 = byu.autospec(ID3,fs,102400)
Gxx4, f4, OASPL4 = byu.autospec(ID4,fs,102400)
Gxx5, f5, OASPL5 = byu.autospec(ID5,fs,102400)

plt.figure()
plt.subplot(311)
plt.plot(f1, Gxx1,'r--')
plt.plot(f2, Gxx2,'b--')
plt.legend(['1st Ambient','2nd Ambient'])
plt.title('Gxx Two Ambient Measurements 11-6-2019, New AFR 11/6 rec. by TC4034')
plt.xlabel('Frequency (Hz)')
plt.grid()
plt.subplot(312)
plt.plot(f3,Gxx3,'b')
plt.title('Gxx 300KHz Sine gen. by TC4038 11/6 rec. by TC4034')
plt.xlabel('Frequency (Hz)')
plt.grid()
plt.subplot(313)
plt.plot(f5,Gxx5,'b')
plt.plot(f4,Gxx4,'r')
plt.legend(['5s Chirp','1s Chirp'])
plt.grid()
plt.title('Gxx Two 100kHz-300KHz Chirps 11/6 rec. by TC4034')
plt.xlabel('Frequency (Hz)')

"""