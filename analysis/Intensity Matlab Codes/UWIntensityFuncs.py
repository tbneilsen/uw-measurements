# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 10:47:13 2021

@author: Corey Dobbs
"""
import numpy as np
import math
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu
import computeBlockFFT as cbFFT


    
    
def TRAD_func(fss, Xss, probe_config, rho=1023, c=1478):
    """
    This function uses the traditional p-p method to calculate
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
        Defaults to water. 1023 kg/m**3
    c : float
        Speed of sound in medium. 
        Defaults to 1478 m/s in water 

    Returns
    -------
    I : Array of float64
        Active time-averaged intensity. 3xN matrix, with each row being the
        values in Cartesian coordinates x, y, and z.
    Q : Array of float64
        Reactive time-averaged intensity. 3xN matrix, with each row being 
        the values in Cartesian coordinates x, y, and z.
    U : Array of float64
        Potential energy density
    T : Array of float64
        Kinetic energy density
    E : Array of float64
        Total energy density (T + U)
    EL : Array of float64
        Lagrangian energy density (T - U)
    I_mag : Array of float64
        Time-averaged active intensity magnitude. Array of scalars.
    I_dir : Array of float64
        Angle (in radians) of the active intensity vector (in x-y plane)
    Q_mag : Array of float64
        Time-averaged reactive intensity magnitude. Array of scalars. 
    Q_dir : Array of float64
        Angle (in radians) of the reactive intensity vector (in x-y plane)
    p0 : Array of complex128
        Complex pressure at center of probe. 
    grad_p : Array of complex128
        Complex pressure gradient. 3xN matrix, one row for each Cartesian
        coordinate. 
    u : Array of complex128
        Particle velocity. 3xN matrix, one row for each Cartesian coordinate.
    z : Array of complex128
        Specific acoustic impedance (p0/u). 3xN vector, one row for each
        Cartesian Coordinate. 
        
    Author: Corey Dobbs
    Transcripted from TRAD_func.m in 
    BYU_Acoustics/Energy-based_Acoustics_Page in
    git.physics.byu.edu
    
    Last modified: 6/25/21
    """
    omega = 2*np.pi*fss
    
    num_CH = np.shape(Xss)[0]
    num_Blocks = np.shape(Xss)[1]
    num_Samples = np.size(Xss)
    
    #The following calculates X and the pseudo inverse of X, used in 
    #calculating a least-squares estimate of grad_p.
    
    #Calculate X and pseudo-inverse of X
    num_pairs = (num_CH**2 - num_CH)//2 #Num of unique microphone pairs
    X = np.zeros((num_pairs, 3))
    
    index = 0
    for i in range(num_CH):
        for j in range(i+1,num_CH):
            X[index,:] = probe_config[j,:] - probe_config[i,:]
            index += 1
    pInvX = np.linalg.pinv(X)

    #Initialize matrices
    I_blocks = np.zeros((num_Blocks,3,num_Samples//2))
    Q_blocks = np.copy(I_blocks)
    usq_blocks = np.copy(I_blocks)
    psq_blocks = np.zeros((num_Blocks, num_Samples//2))
    pdiff = np.zeros((num_pairs, num_Samples//2),dtype = complex)
    
    #The following calculates the TRAD intensities of each block. This is
    #done by calculating grad_p and p0 for each block.
    for i in range(num_Blocks):
        
        p =  np.squeeze(Xss[:,i,:])
        
        #Calculate grad_p
        index = 0
        for j in range(num_CH):
            for k in range(j+1,num_CH):
                pdiff[index,:] = p[k,:] - p[j,:]
                index += 1
              
        
        grad_p = np.matmul(pInvX,pdiff)
           
        #Calculate p0
        centerCH = np.where(probe_config == [0,0,0])
        if len(centerCH) == 1:
            #If there is a microphone at the center of the probe
            p0 = p[centerCH,:] #Complex p0
            P0 = abs(p0) #Abs p0
        else:
            #If there is NOT a microphone at the center of the probe
            p0 = np.mean(p,0) #Complex p0
            P0 = np.mean(abs(p0)) #Abs p0
        
        
        #The function np.tile repeats omega and p0 across x, y and z so they can
        #be multiplied by grad_p element-wise in the intensity calculation.
        omega_vec = np.tile(omega,(3,1))
        p0_vec = np.tile(p0,(3,1))
        
        #Intensity for the current block
        Ic = np.conj(1j/(omega_vec*rho)*grad_p)*p0_vec #Complex intensity
        u = 1j/(omega_vec*rho)*grad_p #Particle velocity
        
        I_blocks[i,:,:] = np.real(Ic) #Real part of intensity
        Q_blocks[i,:,:] = np.imag(Ic) #Imaginary part of intensity
        
        usq_blocks[i,:,:] = abs(u)**2 
        psq_blocks[i,:] = P0**2
    
    #Average intensities across blocks
    I = np.squeeze(np.mean(I_blocks,0))
    Q = np.squeeze(np.mean(Q_blocks,0))
     
    #Calculate energy densities
    psq = np.mean(psq_blocks,0)
    usq = np.squeeze(np.mean(usq_blocks,0))
    U = psq/(2*rho*c**2)
    T_vec = rho/2*usq
    T = np.sum(T_vec,0)
    E = T + U
    EL = T - U  
    
    #Find magnitude of intensity vector
    I_mag = np.squeeze(np.sqrt(np.real(I[0,:])**2 + np.real(I[1,:])**2 + np.real(I[2,:])**2))
    Q_mag = np.squeeze(np.sqrt(np.real(Q[0,:])**2 + np.real(Q[1,:])**2 + np.real(I[2,:])**2))
    
    #Find the angle (in radians) of the intensity vector (in x-y plane)
    Ix = I[0,:]
    Iy = I[1,:]
    I_dir = np.zeros(np.size(Ix))
    Qx = Q[0,:]
    Qy = Q[1,:]
    Q_dir = np.copy(I_dir)
    for i in range(len(Ix)):
        I_dir[i] = math.atan2(Iy[i],Ix[i])
        Q_dir[i] = math.atan2(Qy[i],Qx[i])
    
    z = np.tile(p0,(3,1))/u
    
    return I, Q, U, T, E, EL, I_mag, I_dir, Q_mag, Q_dir, p0, grad_p, u, z

def calcIntensity(fs, L, mic_Dist, filepath, ID_num, c = 1478, rho = 1023, I_ref = 6.61e-19, P_ref = 1e-6):
    """
    This function might be unnecessary. Not very robust. 
    See Intensity_Calculations.py

    Parameters
    ----------
    fs : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    mic_Dist : TYPE
        DESCRIPTION.
    filepath : TYPE
        DESCRIPTION.
    ID_num : TYPE
        DESCRIPTION.
    c : float, optional
        The speed of sound in the used medium. The default is 1478 for water.
    rho : float, optional
        Density of the medium. The default is 1023 for water.
    I_ref : float, optional
        DESCRIPTION. The default is 6.61e-19.
    P_ref : float, optional
        DESCRIPTION. The default is 1e-6.

    Returns
    -------
    I : TYPE
        DESCRIPTION.
    Q : TYPE
        DESCRIPTION.
    U : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    E : TYPE
        DESCRIPTION.
    EL : TYPE
        DESCRIPTION.
    I_mag : TYPE
        DESCRIPTION.
    I_dir : TYPE
        DESCRIPTION.
    Q_mag : TYPE
        DESCRIPTION.
    Q_dir : TYPE
        DESCRIPTION.
    p0 : TYPE
        DESCRIPTION.
    grad_p : TYPE
        DESCRIPTION.
    u : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Author: Corey Dobbs
    Transcripted from the work of Gabriel Fronk in his thesis
    Last Modified: 6/25/21
    """
    
    
    ns =  L*fs #number of samples
    df = fs/ns #sample spacing in frequency domain
    fss = np.transpose(np.arange(0,(fs/2),df)) #single-sided frequency array
    
    #window and weighting functions
    w = np.hanning(ns)
    W = np.mean(w*np.conj(w)) #used to scale the ASD for energy conservation
    
    # Load in data, make sure binfileload.py is downloaded from byuarg library.
    #This is also a good place to account for the amplification factor from 
    #the NEXUS conditioning amplifier (in this case divide by .10 to account 
    #for the 10 mV/Pa amplification).
    x1 = byu.binfileload(filepath,'ID',ID_num,1)/.10
    x2 = byu.binfileload(filepath,'ID',ID_num,2)/.10
    
    Xss1, num_Blocks = cbFFT.computeBlockFFT(x1,ns,w,W)
    Xss2, num_Blocks = cbFFT.computeBlockFFT(x2,ns,w,W)
    
    
    Xss = np.concatenate((np.reshape(Xss1,(1,np.size(Xss1,0),np.size(Xss1,1))),np.reshape(Xss2,(1,np.size(Xss1,0),np.size(Xss1,1))))) 
    
    probe_config = np.array([[0,.12,0],[0,-.12,0]])
        
    I, Q, U, T, E, EL, I_mag, I_dir, Q_mag, Q_dir, p0, grad_p, u, z = TRAD_func(fss,Xss,probe_config)
    
    return fss, I, Q, U, T, E, EL, I_mag, I_dir, Q_mag, Q_dir, p0, grad_p, u, z
    

    
    
    
    