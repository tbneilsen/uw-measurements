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
#from scipy.fftpack import ifft
from datetime import datetime
now = datetime.now()
now_str = now.strftime('%Y_%m_%d')
import os
import pdb
######################################
# Force matplotlib to perform with davids plotting defaults 
# for "paper" or "presentation"

savefolder = 'Out_Figures'
import sys
sys.path.insert(1, "C:\\Users\\cvong\\Box\\UW Research\\Code\\python-plotting-lib")
from fig_defaults import set_defaults as set_fig_defaults
_ = set_fig_defaults("presentation")

######################################
########### Formatting ###############
#loaded_signal = byu.binfileload(path, IDname, IDnum, CHnum)
#number of files if beginning at ID001
#loads as pressures in (Pa)
######################################
files_begin     = 13    #first file (python indexes from 0)
files_end       = 14    #last file (python indexes up to last value)
path = "C:/Users/cvong/Box/UW Research/Underwater Measurements/2020-02-12 Tank Test Sweeps"
fs = 1240000            #sample frequency (Hz)
t_rec_len = 5           #time record length (s)
Ts = 1/fs               #sampling period (s)
N = fs*t_rec_len        #total number of samples
N = int(N) 
n=15                    #number of bins
fmin,fmax = 10000,600000    #for plotting limits
dBmin,dBmax = 40,150        #for plotting limits
# calculate the bin size using n
# note that this is how you can determine how many bins you want to use
# in your averages, you can define your own NS or use this method
ns = int(2 * N / (n + 1))       #samples per block (often use power of 2)
pref = 1e-6                     #reference pressure for water

def prev_power_of_two(x):       #determine power of 2 best for ns
    prev = False
    ii = 0
    while not prev:
        val = 2 ** ii
        if val > x:
            break
        else:
            ii += 1
    return ii - 1
#ns = 2**prev_power_of_two(fs)
                                        #######################################
                                        #set time record window to view...
sampleview = N                          #how many samples you want to focus...  
t_firstwindow = Ts*sampleview           #on."N" being the full value. Adjust...
time = np.arange(0,t_firstwindow,Ts)    #for viewing first window of samples. 
                                        #######################################
print('loading data...')                                     
#load all binfiles for generated (CH0) and received (CH1) signals
all_sig_gen = np.empty((files_end-files_begin,N))
all_sig_rec = np.empty((files_end-files_begin,N))

for ii in range(files_begin,files_end):
    all_sig_gen[ii-files_begin,:] = byu.binfileload(path,"ID",ii,0)   #signals generated
    all_sig_rec[ii-files_begin,:] = byu.binfileload(path,"ID",ii,1)   #signals recorded



print('plotting time waveforms...')
plt.rcParams['agg.path.chunksize']=20000  #needed for plotting long vectors

for ii in range(files_begin,files_end): 
    plt.figure()
    plt.plot(time[:], all_sig_gen[ii-files_begin,:])
    #plt.xlim(0,max(time)/10)
    plt.title(f"T-Waveform for Sig. Gen. {ii}")
    plt.xlabel('Time')
    #pdb.set_trace()
    plt.savefig(os.path.join(savefolder,f't_wave_gen{ii}_{now_str}'))
    
    plt.figure()
    plt.plot(time[:], all_sig_rec[ii-files_begin,:])
    plt.xlim(0.326,0.333) #plt.xlim(0,max(time)/10)
    plt.ylim(-1.5,1.5)
    plt.title(f"T-Waveform for Received Signal {ii}")
    plt.xlabel('Time')
    plt.savefig(os.path.join(savefolder,f't_wave_rec{ii}_{now_str}'))
plt.show()  


"""
plt.axvline(x= ,linestyle = 'dashed',color='r' )    #vertical line to show value
plt.text(xvalue, yvalue, 'text here',color = 'r')
plt.axhline(y= ,linestyle = 'dashed',color='r' )    #horizontal line to show value
plt.text(xvalue, yvalue, 'text here',color = 'r')
     
"""




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

print('data loaded, calc autospec')
Gxx_gen = np.empty((files_end-files_begin,int(ns/2)),dtype = np.complex)
OASPL_gen = np.empty((files_end-files_begin))          #OASPL hard coded for Pref = 20microPa, this is wrong for water
Gxx_rec = np.empty((files_end-files_begin,int(ns/2)),dtype = np.complex)
OASPL_rec = np.empty((files_end-files_begin))
for ii in range(files_begin,files_end):
    Gxx_gen[ii-files_begin], f_new, OASPL_gen[ii-files_begin] = byu.autospec(all_sig_gen[ii-files_begin],fs, ns, N)
    Gxx_rec[ii-files_begin], f_new, OASPL_rec[ii-files_begin] = byu.autospec(all_sig_rec[ii-files_begin],fs, ns, N)

  


for ii in range(files_begin,files_end): 
    plt.figure()
    plt.plot(f_new,Gxx_gen[ii])
#    plt.yscale('log')
    plt.title(f"Autospectra Generated Signal {ii}")
    plt.xlabel('Frequency (Hz)')
    plt.savefig(os.path.join(savefolder,f'Gxx_gen{ii}_{now_str}'))
    
    plt.figure()
    plt.plot(f_new,Gxx_rec[ii])
#    plt.yscale('log')
    plt.title(f"Autospectra Received Signal {ii}")
    plt.xlabel('Frequency (Hz)')
    plt.savefig(os.path.join(savefolder,f'Gxx_rec{ii}_{now_str}'))
plt.show()       


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

"""
print('calculating crossspec')
Gxy = np.empty((files_end-files_begin,int(ns/2)),dtype = np.complex)
for ii in range(files_begin,files_end):
    Gxy[ii-files_begin,:], f_new = byu.crossspec(all_sig_gen[ii-files_begin], all_sig_rec[ii-files_begin], fs, ns, N)

"""

"""    
for ii in range(files_begin,files_end): 
    plt.figure()
    plt.plot(f_new,Gxy[ii])
#    plt.yscale('log')
    plt.title(f"Gxy {ii}")
    plt.xlabel('Frequency (Hz)')
    plt.savefig(os.path.join(savefolder,f'Gxy{ii}_{now_str}'))
plt.show()       
"""  

#########################
#Freq Response Function #
#########################
"""
print('calculating FRF')
H = np.empty((files_end-files_begin,int(ns/2)),dtype = np.complex)
for ii in range(files_begin,files_end):
    H[ii-files_begin,:] = Gxy[ii-files_begin]/Gxx_gen[ii-files_begin]
"""
    
    
"""    
for ii in range(files_begin,files_end): 
    plt.figure()
    plt.plot(f_new,H[ii])
#    plt.yscale('log')
    plt.title(f"FRF {ii}")
    plt.xlabel('Frequency (Hz)')
    plt.savefig(os.path.join(savefolder,f'H{ii}_{now_str}'))
plt.show()       
"""  
    
#########################
# inverse fast fourier  #
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
"""
print('calculating Impulse Response')


# hanning window for determination of scaling factor
#used in autospec to calc fft
ww = np.hanning(ns)
# no need for the conjugate function in the original
# since hanning is always real, we just multiply
W = np.mean(np.multiply(ww, ww))


fft_scale = 1/(2/ns/fs/W)
ImpRes = np.empty((files_end-files_begin,int(ns)))
for ii in range(files_begin,files_end):
    x_conj = np.flip(np.conj(H[ii-files_begin]))
    comb_in = np.append(np.append(H[ii-files_begin],0),x_conj)
    
    comb_in = comb_in[:-1]
    ImpRes[ii-files_begin] = fft_scale * np.real(np.fft.ifft(comb_in/2))
 
"""    
    
"""    
for ii in range(files_begin,files_end): 
    plt.figure()
    plt.plot(ImpRes[ii])
#    plt.yscale('log')
    plt.title(f"Impulse Response {ii}")
    plt.savefig(os.path.join(savefolder,f'ImpRes{ii}_{now_str}'))

plt.show()   
"""
    
    
    



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
    
"""  
print('calculating spectrogram')
ns = int(ns/2**6)

Spec_gen = []               #specgram is 2D and auto builds freq and time array
runOASPL_gen = np.empty((files_end-files_begin))       #not used below because not accurate for underwater
Spec_rec = []
runOASPL_rec = np.empty((files_end-files_begin))

for ii in range(files_begin,files_end):
    Spec_gen_temp,t_gen,f_gen,_ = byu.specgram(all_sig_gen[ii-files_begin],fs,ns)
    Spec_gen.append(Spec_gen_temp)
    Spec_rec_temp,t_rec,f_rec,_ = byu.specgram(all_sig_rec[ii-files_begin],fs,ns)
    Spec_rec.append(Spec_rec_temp)
"""
"""
for ii in range(files_begin,files_end):
    
    plt.figure()
    plt.pcolormesh(t_gen, f_gen, 20*np.log10(np.abs(Spec_gen[ii-files_begin].T)/pref))   #.t to transpose it so time is on bottom
    plt.ylim((fmin,fmax))
    plt.yscale('log')
    plt.title(f"Spectrogram for sig. generated {ii}")
    plt.colorbar()
    plt.clim((dBmin,dBmax))
    plt.savefig(os.path.join(savefolder,f'spec_gen{ii}_{now_str}'))
   
    plt.figure()
    plt.pcolormesh(t_rec, f_rec, 20*np.log10(np.abs(Spec_rec[ii-files_begin].T)/pref))
    plt.ylim((fmin,fmax))
    plt.yscale('log')
    plt.title(f"Spectrogram for sig. recorded {ii}")
    plt.colorbar()
    plt.clim((dBmin,dBmax))
    pdb.set_trace()
    plt.savefig(os.path.join(savefolder,f'spec_rec{ii}_{now_str}'))
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