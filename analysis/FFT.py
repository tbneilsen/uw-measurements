#%%
import numpy as np
from ESAUdata import ESAUdata
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')

from autospec import autospec

path1 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-11/2021-06-11_scan1'
path2 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-11/2021-06-11_scan2'
path3 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-24/2021-06-24_scan7'
desire = [0,1,2,3,4,5,6,7,8,9,10]
channels = [1]
fs = 150000
totalLength = 0.8
N = int(totalLength*fs)
pref = 1e-6
ns = 2**15
unitflag = 0
# %%
_,_,_,_,ch1,_,_ = ESAUdata(path3, desire, channels, N, N)
# %%
num_pos = len(desire)
Gxx = np.zeros((ns//2,num_pos))
_, f, _ = autospec(ch1[:,0], fs, ns, N, unitflag, pref)
oaspl = np.zeros(num_pos)
for i in desire:
    Gxx[:,i], _, oaspl[i] = autospec(ch1[:,i], fs, ns, N, unitflag, pref)
# %%
rel_Gxx_level = np.zeros((ns//2,num_pos - 1))
for i in desire[:num_pos-1]:
    rel_Gxx_level[:,i] = 10*np.log10(Gxx[:,i+1]/Gxx[:,0])
# %%
start = 0 # this was just me finding the index of where 1 kHz is, probably need a better method
stop = 1311 # this was just me finding the index of where 10 kHz is, probably need a better method
#start = 1311 # this was just me finding the index of where 10 kHz is, probably need a better method
#stop = 13107 # this was just me finding the index of where 100 kHz is, probably need a better method
for i in range(len(rel_Gxx_level[0,:])):
    plt.figure()
    plt.title('Frequency vs TL relative to 0.1 m')
    plt.plot(f[start:stop],rel_Gxx_level[start:stop,i], label = 'rel TL ' + str((i+2)/10) + ' m')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Rel TL (dB)')
    plt.grid()
    plt.legend()
    plt.show()
# %%

# %%
