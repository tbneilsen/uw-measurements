# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:15:41 2021

@author: scott
"""

def calcRelativeTransmissionLoss(rec_sig,rec_start,rec_end,fs,signal_duration,leading_and_trailing,step,c):
    '''
    This function calculates the relative Transmission Loss in the Tank.
    Right now it only works for very specific measurements, but hopefully in the future will be a bit more versatile.
    '''
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
        pref = 1e-6 # Pascals
        x = np.square(x)
        OASPL = 10 * np.log10(np.mean(x)/pref**2)
        return OASPL
    # end OASPLcalc

    start_dist = 0.1
    num_scans = len(rec_sig[0,:])

    shiftval = (step/c)*fs
    t = np.linspace(0,(len(rec_sig[:,0]))/fs,len(rec_sig[:,0]))

    dB = np.zeros(num_scans)
    for i in range(num_scans): 
        dB[i] = OASPLcalc(rec_sig[rec_start+int(i*shiftval):rec_end+int(i*shiftval),i])


    rel_TL = dB - dB[0]*np.ones(len(desire))
    Range = [start_dist + step*i for i in range(len(desire))]

    return Range, rel_TL


from ESAUdata import ESAUdata
from readLogFile import readLogFile
import sys
sys.path.append('/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20')

path1 = 'D:/2021-05-18_scan6' # 100 kHz, 10 points
path2 = 'D:/2021-05-18_scan7' # 100 kHz, 10 points
path3 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-20_scan2' # 100 kHz, 100 points
#path3 = 'D:/2021-05-20_scan2' # 100 kHz, 100 points
path4 = 'D:/2021-05-20_scan3' # 50 kHz, 100 points
path5 = 'D:/2021-05-25_scan1' # 100 kHz, 100 points, Ran is on whole time
desire = [i for i in range(100)]
channel = [0,1]
c = 1478
step = 1.0/99.0
rec_start = 500190 # a value chosen by observing the graphs. Need to find a more accurate way of getting this.
rec_end = 500600 # a value chosen by observing the graphs. Need to find a more accurate way of getting this.

_,_,_,fs,signal_duration,_,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path3)
leading_and_trailing = 1.0
N = fs*(leading_and_trailing+signal_duration)
_, _, _, ch0, rec_sig, _, _ = ESAUdata(path3,desire,channel,N)

Range, rel_TL = calcRelativeTransmissionLoss(rec_sig,rec_start,rec_end,fs,signal_duration,leading_and_trailing,step,c)