# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:59:09 2021

@author: Corey Dobbs
"""

import numpy as np
import scipy as sci

def computeBlockFFT(data, ns, w, W):
    """

    Parameters
    ----------
    data : array
        Single vector of the data to which we will apply the block fft.
        Must be in time domain.
    ns : float
        The number of samples in the data
        The block size.
    w : array
        Window function. Ex: Hanning window
        Used to ensure that data starts and ends at same spot,
        so as to appear periodic
    W : float
        The weight. Ex: W = mean(w w*) 
        (average of window array times conjugate window array)

    Returns
    -------
    Xss : array
        The single-sided scaled, fourier-transformed data
    num_Blocks : float
        The number of blocks into which the data is divided.

    Author: Corey Dobbs
    Transcripted from Matlab code computeBlockFFt.m in 
    BYU_Acoustics/Energy-based_Acoustics/energy-based-acoustics
    in git.physics.byu.edu 
    
    Last Modified: 6/25/21
    """
    
    #Make sure data is right shape
    if np.size(data,0) == 1 and np.size(data,1) > 1:
        np.transpose(data)
    
    #Zero mean
    data = data - np.mean(data)
    
    #Get length of data
    N = len(data)
    #Get number of blocks
    num_Blocks = 2*N//ns - 1
    
    #Chunk data
    blockData_1 = np.transpose(np.reshape(data[0:(num_Blocks*ns//2)],(ns//2, num_Blocks)))
    blockData_2 = np.transpose(np.reshape(data[0 + ns//2:(num_Blocks*ns//2 + ns//2)],(ns//2, num_Blocks)))
    
    blockData = np.append(blockData_1, blockData_2)
    
    blockData = blockData*np.tile(w,(num_Blocks,1))
    
    #FFT
    X = sci.fft(blockData)
    
    #Scale data
    Xss = X[:,0:ns//2]*np.sqrt(2/W) #scipy.fft already divides by ns, but matlab does not
    
    return Xss, num_Blocks