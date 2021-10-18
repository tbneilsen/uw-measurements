#%%
from ESAUdata import ESAUdata
from readLogFile import readLogFile
from pathlib import Path

path = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-10-15/2021-10-15_scan2'

desire = [i for i in range(150)]
channels=[0,1]
_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path)

fs = 1e6
measurement_duration = 2.15
N=int(fs*measurement_duration)
_,_,_,tegam,rec_sig,_,_ = ESAUdata(path,desire,channels,N,N)
import matplotlib.pyplot as plt
plt.figure()
plt.plot(tegam[60000:60100,0])
plt.title('TEGAM Monitor')
plt.figure()
plt.plot(rec_sig[60000:60100,0])
plt.title('Received Signal')

# %%
