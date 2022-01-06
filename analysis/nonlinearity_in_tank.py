#%%
from ESAUdata import ESAUdata
from readLogFile import readLogFile
import matplotlib.pyplot as plt
import numpy as np
from TimeGate_UnderwaterTank import gatefunc
from TimeGate_UnderwaterTank import gateValue
from TimeGate_UnderwaterTank import uwsoundspeed
from ESAUpose import ESAUpose
import os

# get data from files
path1 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-06B/2021-08-06_scan1' # 16 points, 100 kHz sine wave, 40 MHz sampling rate
path = path1 # this is the path to be used in the code
numpoints = 16 # change this when you change path!
save_folder = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2022-01'
desire = [i for i in range(numpoints)]
channels=[1]
_,_,temp,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path)
N=int(fs*measurement_duration)
dt = 1/fs

# time array for measurement
t = np.array([dt*i for i in range(N)])
#t = 1000*t # converting to milliseconds instead of seconds
_,_,_,tegam_monitor,rec_sig,_,_ = ESAUdata(path,desire,channels,N,N) # we just want to look at the recorded data and the tegam monitor

# we want to just look at the direct line path so we will first timegate the signal
AEgir_pose,Ran_pose,dd = ESAUpose(path,desire)
c = uwsoundspeed(depth,temp,S=0.03,model='Garrett') # maybe we ought to check that salinity value
tshort = np.zeros(len(desire))
for i in desire:
    tshort[i],_,_,_ = gateValue(AEgir_pose[i], Ran_pose[i], depth, c)
gated_rec_sig = np.zeros((N,len(desire)))
for i in desire:
    gated_rec_sig[:,i]=gatefunc(rec_sig[:,i],fs,tshort[i],leading,tb4=0.1)

# now we can plot the sent signal over the received signal to see if the curve is steepening as sound travels down the tank
start = 0.1 # the beginning of the xlim argument
end = 0.10008000    # the end of the xlim argument
shift = dd/c  # this is the shift required for each position due to delay in sound travel
              # note that dd is an array so shift is also an array. it has units of seconds
'''
for i in desire:
    plt.figure()
    plt.title('Checking for Nonlinearity\n' + 'Source Receiver Distance: ' + str(round(dd[i],3)) + 'm')
    plt.plot(t,rec_sig[:,i])
    plt.xlim(start + shift[i],end + shift[i])
    plt.xlabel('Time (sec)')
    plt.ylabel('Amplitude (volts)')
'''
plt.figure()
plt.title('Checking for Nonlinearity')
for i in [0,15]:
    plt.plot(t,rec_sig[:,i],label = 'SRD: ' + str(round(dd[i],3)) + 'm')
    plt.xlim(start + shift[i],end + shift[i])
plt.xlabel('Time (sec)')
plt.ylim(-1,1)
plt.ylabel('Amplitude (volts)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
save_name = 'nonlinearity_check_time_domain_2points'
plt.savefig(os.path.join(save_folder, save_name))

# %%
