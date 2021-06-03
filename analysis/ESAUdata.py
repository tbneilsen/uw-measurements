# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 19:55:16 2020

This function is to load in the data from scans easily and output the desired 
information needed. With any number of channels desired.  

@author: cvongsaw
"""

def ESAUdata(path,desire=[0],channels=[0],N = 450e3):
    """
    Parameters
    ----------
    path:       string;
                file path name
    desire:     list;
                Desired scan positions to analyze. Should be input as a list of 
                ordered scans. Defaults to only the first measurements data.
    channels:   list, Optional;
                Input channels used (channel numbers for recording data listed)
                Default to [0], being only the first channel on the spectrum cards
                [0,1] would represent the first two channels, [0,2] would represent
                the first and third channel. 
    N:          float, Optional;
                Number of samples. Defaults to 450 kSamples for an fs of 150 kHz 
                and trec of 3 seconds
                
    Returns
    -------
    gen:        float;
                generated signal for each channel outputting signal
    cal:        float;
                calibration measurement for each channel receiving
    ch0:       float;
                all recorded scans desired from ch0
    ch1:       float;
                all recorded scans desired from ch1
    ch2:       float;
                all recorded scans desired from ch2
    ch3:       float;
                all recorded scans desired from ch3 


    Notes
    -----
    Author: Cameron Vongsawad
    
    Does not currently allow for taking in second chassis daisy chained to first.
    
    Now allows for scan without a calibration measurement by checking if a
    calibration file exists before loading cal file. Else cal = zeros and errors
    message about missing cal measurement is printed.
    cal file format updated from cal.bin or cal (1).bin to cal_000.bin. Allows
    for both potential options if referring to older measurements. It also allows
    for nonsequential channel calibration measurements. 
    
    Now allows for non scan measurements and when the file does not contain 
    a generated signal file.
    
    last modified 2/12/2021
    """
    import numpy as np
    import os.path as check
    ##################################################
    #load generated signal and calibration measurement
    ##################################################
    print('loading data...')
    
    isFile_gen = check.isfile(path+ '/signal_out.bin')
    if isFile_gen == True:
        gen = byu.binfileload(path + '/signal_out.bin')
    else:
        gen = np.empty(int(N),dtype = float)
        print('')
        print('Warning: No generated file found. gen = empty')
    
    
    isFile0 = check.isfile(path+ '/cal_000.bin') or check.isfile(path+ '/cal.bin')
    isFile1 = check.isfile(path+ '/cal_001.bin') or check.isfile(path+ '/cal (1).bin')
    isFile2 = check.isfile(path+ '/cal_002.bin') or check.isfile(path+ '/cal (2).bin')
    isFile3 = check.isfile(path+ '/cal_003.bin') or check.isfile(path+ '/cal (3).bin')
    
    isFile = [isFile0,isFile1,isFile2,isFile3]
    
    #load calibration file for each channel recorded into 1 of 4 columns in the 
    #array for the 4 channels allowed with the spectrum cards. 
    cal = np.empty((len(gen),len(isFile)),dtype = float)
    #for idx,ch in enumerate(channels):
    for idx in range(len(isFile)):
        if isFile[idx] == True:
            if check.isfile(path+ f'/cal_00{idx}.bin') == True:
                cal0 = byu.binfileload(path + f'/cal_00{idx}.bin')
                cal[:,idx] = cal0
            elif check.isfile(path+ f'/cal ({idx}).bin'):
                cal0 = byu.binfileload(path + f'/cal ({idx}).bin')
                cal[:,idx] = cal0
            elif check.isfile(path+ '/cal.bin') == True: 
                cal0 = byu.binfileload(path + '/cal.bin')
                cal[:,idx] = cal0
        else: 
            #cal[:,ch] = np.zeros(len(gen),dtype = float)
            print('')
            print(f'Warning: no ch{idx} calibration file found')
    
    if (cal ==  np.empty((len(gen),4),dtype = float)) is True:
        print('')
        print('Warning: recording error: calibration not recorded, file is empty')
                
                
    """ 
    This is the old version of the cal file loading code. This did not allow for 
    measurements to be taken on unconsecutive channels or without first using ch0
    if isFile0 == True:
        if check.isfile(path+ '/cal_000.bin') == True:
            cal0 = byu.binfileload(path + '/cal_000.bin')
            cal[:,0] = cal0
        else: 
            cal0 = byu.binfileload(path + '/cal.bin')
            cal[:,0] = cal0
    else: 
        cal0 = np.zeros(len(gen),dtype = float)
        print('')
        print('Warning: no ch0 calibration file found')
    
    if (cal0 == np.zeros(len(gen),dtype = float)) is True:
        print('')
        print('Warning: recording error: calibration not recorded, file is empty')
        
    if isFile1 == True:
        if check.isfile(path+ '/cal_001.bin') == True:
            cal1 = byu.binfileload(path + '/cal_001.bin')
            cal[:,1] = cal1
        else: 
            cal1 = byu.binfileload(path + '/cal (1).bin')
            cal[:,1] = cal1
    if isFile2 == True:
        if check.isfile(path+ '/cal_002.bin') == True:
            cal2 = byu.binfileload(path + '/cal_002.bin')
            cal[:,2] = cal2
        else: 
            cal2 = byu.binfileload(path + '/cal (2).bin')
            cal[:,2] = cal2
    if isFile3 == True:
        if check.isfile(path+ '/cal_003.bin') == True:
            cal3 = byu.binfileload(path + '/cal_003.bin')
            cal[:,3] = cal3
        else: 
            cal3 = byu.binfileload(path + '/cal (3).bin')
            cal[:,3] = cal3
            """    
    #load all scan binfiles
    ch0 = np.empty((len(gen),len(desire)))
    ch1 = np.empty((len(gen),len(desire)))
    ch2 = np.empty((len(gen),len(desire)))
    ch3 = np.empty((len(gen),len(desire)))
    #pull only the "desire" values and their index value from the list to populate array
    #to view channel 0
    if 0 in channels:
        for idx,ich in enumerate(desire):
            ch0[:,idx] = byu.binfileload(path,"ID",ich,0)
    #to view channel 1
    if 1 in channels:
        for idx,ich in enumerate(desire):
            ch1[:,idx] = byu.binfileload(path,"ID",ich,1)
    #to view channel 2
    if 2 in channels:
        for idx,ich in enumerate(desire):
            ch2[:,idx] = byu.binfileload(path,"ID",ich,2)
    #to view channel 3
    if 3 in channels:
        for idx,ich in enumerate(desire):
            ch3[:,idx] = byu.binfileload(path,"ID",ich,3)   
        
    return gen, cal, ch0, ch1, ch2, ch3





"""
https://www.computerhope.com/issues/ch001721.htm

Create Function to Read Log Files and find the values of:
    -Cal Pose
    -Specific Scan pose
    -fs
    -number of sampels
    -sigL (leading zeros, signal, trailing zeros)
    -Freq Band of Signal & Type of signal
    -Temperature
    -Water Depth
"""