# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 10:47:13 2021

@author: Corey Dobbs
"""
import numpy as np
import scipy as sci
import sys
sys.path.append("path to byu arg lib")
import byuarglib as byu

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
    w : float? array?
        Window
    W : TYPE
        The weight. Ex: W = mean(W**2)

    Returns
    -------
    Xss : array
        The single-sided scaled data
    num_Blocks : float
        The number of blocks into which the data was divided.

    Author: Corey Dobbs
    Transcripted from Matlab code computeBlockFFt.m in 
    BYU_Acoustics/Energy-based_Acoustics/energy-based-acoustics
    in git.physics.byu.edu 
    
    Last Modified: 6/15/21
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
    Xss = X[:,1:np.floor(ns/2)]/ns*np.sqrt(2/W)
    
    return Xss, num_Blocks 
    
    
def TRAD_func(fss, Xss, probe_config, rho=997, c=1478):
    """
    This function uses the traditional method to calculate
    estimated intensity vectors from multi-microphone intensity
    probes given single-sided fft's.
    
    TRAD - short for traditional method

    Parameters
    ----------
    fss : array
        Single-sided frequency array corresponding to Xss
    Xss : array
        matrix of blocked, single-sided fft's.
        Organized as [CH_number, block_number, sample]
        Ex: fft of the 2nd microphone's 15th block would be
        Xss[2,15,:]
    probe_config : array
        An nx3 matrix, where n is the number of of channels
        per probe. It contains the location of each microphone 
        relative to the geometric center of the probe, in 
        Cartesian coordinates. If there is a center microphone,
        it would have the coordinates (0,0,0)
    rho : float
        Density of medium
        Defaults to water. 997 kg/m**3
    c : float
        Speed of sound in medium. 
        Defaults to 1478 m/s in water 

    Returns
    -------
    I : TYPE
        Active time-averaged intensity
    Q : TYPE
        Reactive time-averaged intensity
    U : TYPE
        Potential energy density (same for TRAD and PAGE)
    T : TYPE
        Kinetic energy density
    E : TYPE
        Total energy density (T + U)
    EL : TYPE
        Lagrangian energy density (T - U)

    Author: Corey Dobbs
    Transcripted from TRAD_func.m in 
    BYU_Acoustics/Energy-based_Acoustics_Page in
    git.physics.byu.edu
    
    Last modified: 6/15/21
    """
    
    omega = 2*np.pi*fss
    
    size_fft = np.shape(Xss)
    num_CH = size_fft[0]
    num_Blocks = size_fft[1]
    num_Samples = size_fft[2]
    
    
    #The following calculates X and the pseudo inverse of X, used in 
    #calculating a least-squares estimate of grad_p.
    
    #Calculate X and pseudo-inverse of X
    num_pairs = (num_CH**2 - num_CH)/2 #Num of unique microphone pairs
    X = np.zeros(num_pairs, 3)
    
    index = 0
    for i in range(num_CH):
        for j in range(i+1,num_CH):
            X[index,:] = probe_config[j,:] - probe_config[i,:]
            index += 1
    pInvX = np.linalg.pinv(X)

    #Initialize matrices
    I_blocks = np.zeros((num_Blocks,3,num_Samples/2))
    Q_blocks = np.copy(I_blocks)
    usq_blocks = np.copy(I_blocks)
    psq_blocks = np.zeros((num_Blocks, num_Samples/2))
    pdiff = np.zeros((num_pairs, num_Samples/2))
    
    #The following calculates the TRAD intensities of each block. This is
    #done by calculating grad_p and p0 for each block.
    for i in range(num_Blocks):
        p =  np.squeeze(Xss[:,i,:])
        
        #Calculate grad_p
        index = 0
        for j in range(num_CH):
            for k in range(j+1,num_CH):
                pdiff[index,:] = p[k,:]
                index += 0
                
        grad_p = pInvX*pdiff
            
        #Calculate p0
        centerCH = probe_config.index(0)
    return I, Q, U, T, E, EL


    