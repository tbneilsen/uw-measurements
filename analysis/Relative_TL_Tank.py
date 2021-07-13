# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:15:41 2021

@author: scott
"""
#%%

def calcRelativeTransmissionLoss(rec_sig,rec_start,rec_end):
    '''
    This function calculates the relative Transmission Loss in the Tank.
    Right now it only works for very specific measurements, but hopefully in the future will be a bit more versatile.

    Parameters
    ----------
    rec_sig : ndarray of float
        The recorded signal that you want to find the relative transmission loss for. This should be a 2 dimmensional
        array of shape (N,num_pos) where N is the number of samples and num_pos is how many positions were in your 
        scan. If num_pos = 1 then you can't find relative transmission loss and this function is not for you.
    rec_start : ndarray of int
        These are the indices 
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
        pref = 12.34684596708469078406978 # Pascals
        #x = np.square(x)
        #OASPL = 10 * np.log10(np.mean(x)/pref**2)
        OASPL = 10 * np.log10(np.var(x)/pref**2)
        return OASPL
    # end OASPLcalc

    dB = np.zeros(len(rec_sig[0,:]))

    for i in range(len(rec_sig[0,:])): 
        dB[i] = OASPLcalc(rec_sig[rec_start[i]:rec_end[i],i])
    
    rel_TL = dB - dB[0]*np.ones(len(desire))

    return rel_TL


from ESAUdata import ESAUdata
from ESAUpose import ESAUpose
from readLogFile import readLogFile
import TimeGate_UnderwaterTank as TG
import numpy as np
import sys
sys.path.append('/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20')

path3 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-20_scan2' # 100 kHz, 100 points
path4 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-20_scan3' # 50 kHz, 100 points
path5 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-25_scan1' # 100 kHz, 100 points, Ran is on whole time
path6 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-07-09/2021-07-09_scan' # 100 kHz, 100 points, signal on longer
path7 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-07-09/2021-07-09_scan1' # 100 kHz, 100 points, signal on longer, longer settling time
path_used = path7
desire = [i for i in range(100)]
channel = [0,1]
c = 1478


_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path_used)
A, R, dd = ESAUpose(path_used,desire)
tshort = np.zeros(len(desire))
tside = np.zeros(len(desire))
tdirect = np.zeros(len(desire))
for i in desire:
    tshort[i],tside[i],tdirect[i],_ = TG.gateValue(A[i], R[i], depth, c, 'tank',False)

N = int(fs*measurement_duration)

_, _, _, _, rec_sig, _, _ = ESAUdata(path_used,desire,channel,N)

gated_rec_sig = np.zeros((N,len(desire)))
for i in desire:
    gated_rec_sig[:,i] = TG.gatefunc(rec_sig[:,i],fs,leading + tside[i], 0.1)

import matplotlib.pyplot as plt

rec_start = fs*(tdirect + leading)
rec_end = fs*(tside + leading)

rec_start = rec_start.astype('int')
rec_end = rec_end.astype('int')

t = np.linspace(0,(N-1)/fs,N)


'''
for i in desire:
    plt.figure()
    plt.plot(t[rec_start[i]:rec_end[i]],gated_rec_sig[rec_start[i]:rec_end[i],i], label = str(np.round(dd[i],4)) + ' m')
    plt.ylim(-2.5,2.5)
    plt.legend()
    plt.show()
'''

rel_TL = calcRelativeTransmissionLoss(gated_rec_sig,rec_start,rec_end)
plt.figure()
plt.plot(dd,rel_TL)
plt.show()

# %%
