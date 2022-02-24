# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:59:09 2021

@author: Corey Dobbs
"""

import numpy as np
import numpy.matlib
import scipy as sci

def computeBlockFFT(data, ns):
    """

    Parameters
    ----------
    data : array
        Single vector of the data to which we will apply the block fft.
        Must be in time domain.
    ns : int
        The block size, i.e. number of samples in each block.
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
        The single-sided scaled, Fourier-transformed data
    num_Blocks : float
        The number of blocks into which the data is divided.

    Author: Corey Dobbs
    Transcripted from Matlab code computeBlockFFt.m in 
    BYU_Acoustics/Energy-based_Acoustics/energy-based-acoustics
    in git.physics.byu.edu 
    
    Last Modified: 11/9/21
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
    
    #Make sure window function is same size as ns
    w = np.hanning(ns)
    W = np.mean(w*np.conj(w))
    
    #Chunk data
    
    #get the array indices where each block is located
    blockmat = (np.matlib.repmat(np.arange(1, ns+1), num_Blocks, 1)
                + np.matlib.repmat(
        ((ns/2) * np.arange(num_Blocks)[:, np.newaxis]).astype(np.int),
        1, ns) - 1).astype(np.uint32)

    
    #Make array of size (num_Blocks, ns)
    #build the matrix, and multiply by the window function
    blockData = np.multiply(np.squeeze(data[blockmat]), w)
        
    #FFT
    X = sci.fft(blockData)
    
    #Scale data
    Xss = X[:,0:ns//2]*np.sqrt(2/W) #scipy.fft already divides by ns, but matlab does not
    
    return Xss, num_Blocks