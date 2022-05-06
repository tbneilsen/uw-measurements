import os
import pickle
filename_bellhop = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/bellhop/defaultoutput/2022-02/modeling_paper_bellhop_pickle'
infile = open(filename_bellhop,'rb')
bellhop_dict = pickle.load(infile)
infile.close()
bellhop_src_depth = bellhop_dict['src_depth']
bellhop_rec_depth = bellhop_dict['rec_depth']
bellhop_ranges = bellhop_dict['ranges']
bellhop_frequency = bellhop_dict['frequency']
bellhop_tl = bellhop_dict['transmission loss']
import matplotlib.pyplot as plt
import numpy as np

def findIndex(array,val):
  # This function takes an array and a value. This value should be in the array somewhere, you just don't know where
  # It returns the index of the closest thing to that value
  # WARNING: if you put a value in that is not inside that array, it will still return the closest thing to it
  # EXAMPLE: array = [6,7,8,9,10,11,12], val = 7, function returns: 1
  # EXAMPLE: array = [6.002,7.123,8.34,9.01,10.134,11.12,12.00101], val = 7, function returns: 1
  # BAD EXAMPLE: array = [6,7,8,9,10,11,12], val = 18000, function returns: 6
    import numpy as np
    adjust_array = np.array(array) # turns into a numpy array. A bit easier to work with
    adjust_array = adjust_array - val # This is subtracting your value from each element in the array
    adjust_array = abs(adjust_array) # makes all values in adjust array positive
    minimum = min(adjust_array) # finds the smallest value in the array (should be pretty close to 0 if not 0)
    index = np.where(adjust_array == minimum)[0][0] # finds the index of that smallest value
    return index

def relOAPSLfft(Gxx,freq_index,plusorminus):
    listlen = len(Gxx[0,:])
    OASPL = np.zeros(listlen)
    for i in range(listlen):
        OASPL[i] = 10*np.log10(np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,i])/np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,0]))
    return OASPL

from ESAUdata import ESAUdata
from readLogFile import readLogFile
from ESAUpose import ESAUpose

import sys
#sys.path.append('/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20')
sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')

from autospec import autospec

path15 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-10-15/2021-10-15_scan2' # 100 kHz 150 points steady state
path17 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan4' # 71 kHz 150 points long sine wave
path18 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan5' # 82 kHz 150 points long
path19 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan8' # 115 kHz 150 points long
path20 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan9' # 133 kHz 150 points long
path21 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan10' # 200 kHz 150 points long

path_used = path21
f = 200000.0

desire = [i for i in range(150)]
channel = [0,1]
c =1486

_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path_used)
A, R, dd = ESAUpose(path_used,desire)
N = int(fs*measurement_duration)
_, _, _, _, rec_sig, _, _ = ESAUdata(path_used,desire,channel,N)

ns = 2**15
unitflag = 0
pref = 1e-6
Gxx = np.zeros((ns//2,len(desire)))
_, freqs, _ = autospec(rec_sig[:,0], fs, ns, N, unitflag, pref)
for i in desire:
    Gxx[:,i], _, _ = autospec(rec_sig[:,i], fs, ns, N, unitflag, pref)

from Relative_TL_Tank import calcRelativeTransmissionLoss

plusorminus = 50
bellhop_index = findIndex(bellhop_frequency,f)
freq_index = findIndex(freqs,f)
rel_TL_fft = relOAPSLfft(Gxx,freq_index,plusorminus)

SAVE_FOLDER = "/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2022-02"

plt.figure()
# You also need to comment/uncomment out plt.savefig() in or out of loop!!
plt.plot(bellhop_ranges, bellhop_tl[0,0,:,bellhop_index] - bellhop_tl[0,0,0,bellhop_index]*np.ones(len(bellhop_tl[0,0,:,bellhop_index])),\
        label = 'Bellhop') #"calc TL re " + str(round(ranges[0],3)) + ' m @ ' + str(f) + " Hz") # Calc Relative TL
# this is where you can plot actual data over the calculated tl
plt.plot(dd,rel_TL_fft, label = 'measured') # true Relative TL
plt.xlabel('Range, m')
plt.ylabel("TL, dB re "+str(round(dd[0],3))+" m") # change re depending on what it's relative to
plt.title('Range vs Relative Transmission Loss\
    \nReciever Depth: ' + str(bellhop_rec_depth[0]) + ' m' +\
    '\nSource Depth: ' + str(bellhop_src_depth[0]) + ' m' +\
    '\nFrequency: ' + str(bellhop_frequency[bellhop_index]) + ' Hz')
# plt.title('Short Sine Wave Time Gated \n '+str(f)+' Hz')
plt.grid()
plt.legend()
plt.gcf().tight_layout()
save_name = 'comp_bellhop_' + 'range_vs_rtl_' + '@freq' + str(f) + 'Hz' + '.png'
#plt.savefig(os.path.join( SAVE_FOLDER, save_name))
plt.savefig(os.path.join( SAVE_FOLDER, save_name))