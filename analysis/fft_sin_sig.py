#%%
import sys
import numpy as np

def findFreqIndex(f_array,freq): # function that finds the closest index to a desired frequency
        adjust_array = abs(f_array - freq)
        val = min(adjust_array)
        index = np.where(adjust_array == val)[0][0]
        return index

def relOAPSLfft(Gxx,freq_index,plusorminus):
    listlen = len(Gxx[0,:])
    OASPL = np.zeros(listlen)
    for i in range(listlen):
        OASPL[i] = 10*np.log10(np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,i])/np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,0]))
    return OASPL

def relOASPLtime(rec_sig):
    listlen = len(rec_sig[0,:])
    OASPL = np.zeros(listlen)
    for i in range(listlen):
        OASPL[i] = 10*np.log10(np.var(rec_sig[:,i])/np.var(rec_sig[:,0]))
    return OASPL

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')

from autospec import autospec
from ESAUdata import ESAUdata
from readLogFile import readLogFile

path10 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan' # 100 kHz
path11 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan1' # 80 kHz
path12 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan2' # 50 kHz
path13 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan3' # 30 kHz
path14 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan4' # 10 kHz

path_used = path14
freq = 10000
desire = [i for i in range(100)]
channel = [0,1]

_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path_used)
N = int(fs*measurement_duration)
_, _, _, _, rec_sig, _, _ = ESAUdata(path_used,desire,channel,N)
ns = 2**15
unitflag = 0
pref = 1e-6
Gxx = np.zeros((ns//2,len(desire)))
OASPL = np.zeros(len(desire))
_, f, _ = autospec(rec_sig, fs, ns, N, unitflag, pref)
for i in desire:
    Gxx[:,i], _, OASPL[i] = autospec(rec_sig[:,i]-np.mean(rec_sig[:,i]), fs, ns, N, unitflag, pref)

import matplotlib.pyplot as plt

index = findFreqIndex(f,freq)
plusorminus = 50

'''
for i in [0,10,20,50,99]:
    plt.figure()
    plt.plot(f[index-plusorminus:index+plusorminus],Gxx[index-plusorminus:index+plusorminus,i])
    #plt.ylim((0,5e-7))
    plt.title(str(freq) + ' Hz\n'\
        +'Range: ' + str(round(0.1+i/99,3)) + ' m\n'\
         + 'Plus or minus ' + str(plusorminus) + ' bins')
    plt.ylabel('Gxx')
    plt.xlabel('frequency (Hz)')
    plt.show()
'''
rel_OASPL_fft = relOAPSLfft(Gxx,index,plusorminus)
rel_OASPL_time = relOASPLtime(rec_sig)
dd = np.array([0.1 + round(i/99,3) for i in range(100)])
plt.figure()
plt.plot(dd,rel_OASPL_fft, label = 'from fft')
plt.plot(dd,rel_OASPL_time, label = 'from time domain')
plt.legend()
plt.title('Relative Transmission Loss\n'\
    + 'Relative to 0.1 meters\n'\
        + str(freq) + ' Hz\n')
plt.ylabel('dB re 0.1 meters')
plt.xlabel('range (meters)')
plt.ylim((-35,6))
plt.grid()
plt.show()
# %%
