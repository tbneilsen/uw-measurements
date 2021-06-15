# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 10:47:13 2021

@author: Corey Dobbs
"""
import numpy as np
import scipy as sci

def computeBlockFFT(data, ns, w, W):
    """

    Parameters
    ----------
    data : array
        single vector of the data to which we will apply the block fft
    ns : float
        The number of samples in the data
    w : float? array?
        Window
    W : TYPE
        DESCRIPTION.

    Returns
    -------
    Xss : TYPE
        DESCRIPTION.
    num_Blocks : float
        DESCRIPTION.

    """
    
    #Make sure data is right shape
    if np.size(data,0) == 1 and np.size(data,1) > 1:
        np.transpose(data)
    
    #Zero mean
    data = data - np.mean(data)
    
    #Get length of data
    N = len(data)
    #Get number of blocks
    num_Blocks = np.floor(2*N/ns - 1)
    #Chunk data
    blockData_1 = np.transpose(np.reshape(data[0:(num_Blocks*ns/2)],[ns/2, num_Blocks]))
    blockData_2 = np.transpose(np.reshape(data[0 + ns/2:(num_Blocks*ns/2 + ns/2)],[ns/2, num_Blocks]))
    
    blockData = [blockData_1, blockData_2]
    blockData = blockData*np.tile(w,(num_Blocks,1))
    
    #FFT
    X = sci.fft(blockData)
    
    #Scale data
    Xss = np.copy(X)
    
    return Xss, num_Blocks 
    
    
def TRAD_func():
    
    x = 1
    return x


    