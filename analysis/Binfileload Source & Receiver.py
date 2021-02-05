# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:28:12 2020

BinFileLoad Source and Receiver

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
path = "C:/Users/cvong/Box/UW Research/Underwater Measurements/2020-02-13 Tank Test Sweeps 2"
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