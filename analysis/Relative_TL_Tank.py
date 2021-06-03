# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:15:41 2021

@author: scott
"""

import numpy as np
def OASPLcalc(x):
    """
    Calculate the overall sound pressure level of a signal
    This function takes a pressure waveform, x, and calculates the overall
    sound pressure level.
    Parameters
    ----------
    x : ndarray of float
        Time series pressure signal. Should be a single dimensional array.
    Returns
    -------
    OASPL : float
        OASPL is the calculated overall sound pressure level
    Notes
    -------
    Author: David Van Komen (david.vankomen@gmail.com)
    Last Modified: 07/02/18
    Based on: OASPLcalc.m by Kent Gee
    """
    # square the pressure array
    x = np.square(x)
    OASPL = 10 * np.log10(np.mean(x)/2e-5**2)
    return OASPL
# end OASPLcalc

from ESAUdata import ESAUdata
from readLogFile import readLogFile
from TimeGate_UnderwaterTank import gateValue
from TimeGate_UnderwaterTank import gatefunc

path1 = 'D:/2021-05-18_scan6' # 100 kHz, 10 points
path2 = 'D:/2021-05-18_scan7' # 100 kHz, 10 points
path3 = 'D:/2021-05-20_scan2' # 100 kHz, 100 points
path4 = 'D:/2021-05-20_scan3' # 50 kHz, 100 points
path5 = 'D:/2021-05-25_scan1' # 100 kHz, 100 points, Ran is on whole time
desire = [i for i in range(100)]
channel = [0,1]
c = 1478

xSource = np.zeros(len(desire))
ySource = np.zeros(len(desire))
zSource = np.zeros(len(desire))
xRec = np.zeros(len(desire))
yRec = np.zeros(len(desire))
zRec = np.zeros(len(desire))
tside = np.zeros(len(desire))
freqMin,freqMax,tempWater,fs,signalDuration,hWater,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path3)
for i in range(10):
    for j in range(10):
        _,_,_,_,_,_,xSource[i*10+j],ySource[i*10+j],zSource[i*10+j],xRec[i*10+j],yRec[i*10+j],zRec[i*10+j] = readLogFile(f'/ID0{i}{j}_001log.txt',path3)
        _,tside[i*10+j],_ = gateValue((xSource[i*10+j],ySource[i*10+j],zSource[i*10+j]), (xRec[i*10+j],yRec[i*10+j],zRec[i*10+j]), hWater, c, 'tank')
# for location just use path1, path2, etc.

N = fs*(1+signalDuration)
leading_zeros = 0.500093
sig_start = int(leading_zeros*fs)
sig_end = sig_start+int(signalDuration*fs)+1
rec_start = 500190
rec_end = 500600
start_dist = 0.1

_, _, _, ch0, ch1, _, _ = ESAUdata(path3,desire,channel,N)

ch1gated = np.zeros((int(N),len(desire)))
tside = tside + leading_zeros

for i in desire:
    ch1gated[:,i] = gatefunc(ch1[:,i],fs,tside[i],0.1)
    
import matplotlib.pyplot as plt

step = 1.0/99.0
shiftval = (step/c)*fs
t = np.linspace(0,(len(ch0[:,0]))/fs,len(ch0[:,0]))

plt.figure()
plt.plot(t,ch0[:,i])
plt.title("Output Signal")
plt.show()

plt.figure()
plt.plot(t[sig_start:sig_end],ch0[sig_start:sig_end,0])
plt.title("Output Signal Zoomed in")
plt.show()

for i in desire:
    plt.figure()
    plt.plot(t[rec_start+int(i*shiftval):rec_end+int(i*shiftval)],ch1[rec_start+int(i*shiftval):rec_end+int(i*shiftval),i])
    plt.title(f"Received Signal Zoomed in (pos {i+1})")
    plt.ylim(-3,3)
    plt.show()

dB = np.zeros(len(desire))
for i in desire:
    dB[i] = OASPLcalc(ch1[rec_start+int(i*shiftval):rec_end+int(i*shiftval),i])


relTL = dB - dB[0]*np.ones(len(desire))
Range = [start_dist + step*i for i in range(len(desire))]

plt.figure()
plt.plot(Range,relTL)
plt.grid()
plt.xlabel("Range in meters")
plt.ylabel("Transmission Loss in dB")
plt.ylim(-25,0)
plt.title("Relative TL vs. Range (using 0.1 m as reference)")
plt.show
