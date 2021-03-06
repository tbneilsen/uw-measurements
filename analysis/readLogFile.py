# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:11:18 2021

@author: Corey Dobbs
"""


def readLogFile(filename,location):
    """
    Code Description:
        This code takes the log file output by ESAU and extracts important details.
        Details described in Returns section.
        The bandwidth, source and receiver positions, sampling frequency, water temperature,
        and water depth are all printed at the end of the code.
        Relevant variables are listed as returns.


    Parameters
    ----------
    filename : String
        The name of the file (type .txt) from which you will pull the experiment parameters.
        Start and end with quotes. 
        Ex: 'ID001_001log.txt'
        
    location : String
        Directory of file that you want to pull
        Ex: 'D:/uw-acoustics-research/uw-meas-codes/underwater-measurements/analysis/'
        IMPORTANT: Check direction of slashes 
        
    Returns
    -------
    freqMin : Float
        The lowest frequency on the frequency range of the sweep. (Hz)
    freqMax : Float
        The highest frequency on the frequency range of the sweep. (Hz)
    tempWater : Float
        Temperature of the water at the time of measurement (degrees Celsius)
    fs : Float
        Sampling frequency of measurement (Hz)
    leadingZeros : Float
        leading zeros before the measurement starts (s)
    signalDuration : Float
        Length of signal (s)
    trailingZeros : Float
        trailing zeros after the signal has been sent (s)
    measurementDuration : Float
        Length of Measurement (s)
    hWater : Float
        Height/depth of water in tank at time of measurement, measured in meters
        from the bottom of the tank
    xSource : Float
        X-position of source (m)
    ySource : Float
        Y-position of source (m)
    zSource : Float
        Z-position of source (m)
    xRec : Float
        X-position of receiver (m)
    yRec : Float
        Y-position of receiver (m)
    zRec : Float
        Z-position of receiver (m)
        
    Notes
    -------
    Created by Corey Dobbs
    
    Last updated 4/8/2021
    """
    
    #Call desired directory
    #import sys
    #sys.path.insert(1,location)
    
    mylines = []
    with open(location+filename,"rt") as myfile: 
        for myline in myfile:
            mylines.append(myline.rstrip('\n'))


    #Find receiver position
    substrR = "Receiver: "
    for line in mylines:          # string to be searched
        index = 0                   # current index: character being compared
        prev = 0                    # previous index: last character compared
        while index < len(line):    # While index has not exceeded string length,
            index = line.find(substrR, index)  # set index to first occurrence of substring
            if index == -1:           # If nothing was found,
                break                   # exit the while loop. 
            receiverPos = "(" + line[index+len(substrR):index+len(substrR)+17] + ")"
                                                 
            prev = index + len(substrR)       # remember this position for next loop.
            index += len(substrR)      # increment the index by the length of substr.
                                      # (Repeat until index > line length)
                              
    #Find source position
    substrS = "Source: "                  
    for line in mylines:          
        index = 0                 
        prev = 0                   
        while index < len(line):   
            index = line.find(substrS, index) 
            if index == -1:         
                break                 
            sourcePos = "(" + line[index+len(substrS):index+len(substrS)+17] + ")"
                                             
            prev = index + len(substrS)      
            index += len(substrS)      
                             
                              
    #Find water depth
    substrD = "Level: "                  
    for line in mylines:         
        index = 0                 
        prev = 0                    
        while index < len(line):   
            index = line.find(substrD, index) 
            if index == -1:          
                break                  
            waterDepth = line[index+len(substrD):index+len(substrD)+10]
                                             
            prev = index + len(substrD)     
            index += len(substrD)
                           
                              
    #Find sampling frequency
    substrF = "Frequency: "                  
    for line in mylines:         
        index = 0                 
        prev = 0                
        while index < len(line):   
            index = line.find(substrF, index) 
            endIndex = line.find("(Hz)")
            if index == -1:           
                break                   
            samplingFreq = line[index+len(substrF):endIndex] + "Hz"
                                             
            prev = index + len(substrF)      
            index += len(substrF)     
                             
                              
    #Find bandwidth
    substrBmin = "from " 
    substrBmax = "to "                 
    for line in mylines:
        index= 0
        prev = 0   
    
        if substrBmin in line:
               a = line.find(substrBmin, index)
               endIndex = line.find('.00')
               fmin = line[a+len(substrBmin):endIndex]
               nextIndex = line.find(substrBmax,index)
               endIndex2 = line.index('.00',nextIndex)
               fmax = line[nextIndex+3:endIndex2]
               break
               if index == -1:           
                break                  
               index += 1
            
               fmin = line[index+len(substrBmin):endIndex]                  
               
    bandwidth = fmin + "-" + fmax + " Hz"


                  

    #Find water temperature
    substrT = "Temp: "                  
    for line in mylines:          
        index = 0                   
        prev = 0                  
        while index < len(line):    
            index = line.find(substrT, index) 
            if index == -1:         
                break                  
            waterTemp = line[index+len(substrT):index+len(substrT)+8]
                                             
            prev = index + len(substrT)       
            index += len(substrT)     
                             
    #Find Leading 0's
    substrLead = "Leading 0's: "
    for line in mylines:
        index = 0
        prev = 0
        while index < len(line):
            index = line.find(substrLead, index)
            endIndex = line.find("  Signal length:")
            if index == -1:    # If nothing was found,
                break           # exit the while loop.
            leading = line[index+len(substrLead):endIndex]

            prev = index + len(substrLead)      # remember this position for next loop.
            index += len(substrLead)        # increment the index by the length of substr.
                                 # (Repeat until index > line length)

    #Find signal length
    substrL = "length: "                  
    for line in mylines:          
        index = 0                  
        prev = 0                    
        while index < len(line):   
            index = line.find(substrL, index) 
            endIndex = line.find("   Trailing")
            if index == -1:           # If nothing was found,
                break                   # exit the while loop. 
            signalLength = line[index+len(substrL):endIndex]
                                             
            prev = index + len(substrL)       # remember this position for next loop.
            index += len(substrL)      # increment the index by the length of substr.
                              # (Repeat until index > line length)

    #Find Trailing 0's
    substrTrail = "Trailing 0's: "
    for line in mylines:
        index = 0
        prev = 0
        while index < len(line):
            index = line.find(substrTrail, index)
            if index == -1:    # If nothing was found,
                break           # exit the while loop.
            trailing = line[index+len(substrTrail):] # no endIndex, needs to go to the end of line

            prev = index + len(substrTrail)      # remember this position for next loop.
            index += len(substrTrail)        # increment the index by the length of substr.
                                 # (Repeat until index > line length)
                              
    print("Bandwidth: ",bandwidth,'\n',"Water Temp: ",waterTemp,'\n',"Source Position: ",sourcePos,'\n',"Receiver Position: ",receiverPos,'\n',"Sampling Freq: ",samplingFreq,'\n',"Water Height: ",waterDepth)
    
    #Convert strings to variables (floats)
    freqMin = float(fmin)
    freqMax = float(fmax)
    tempWater = float(waterTemp[0:-2])
    fs = float(samplingFreq[0:-2])
    Ns = float(signalLength)
    leadingZeros = float(leading)/fs
    signalDuration = Ns/fs
    trailingZeros = float(trailing)/fs
    measurementDuration = leadingZeros + signalDuration + trailingZeros
    hWater = float(waterDepth[0:-1])
    xSource = float(sourcePos[1:6])
    ySource = float(sourcePos[7:12])
    zSource = float(sourcePos[13:-1])
    xRec = float(receiverPos[1:6])
    yRec = float(receiverPos[7:12])
    zRec = float(receiverPos[13:-1])
    
    return freqMin,freqMax,tempWater,fs,leadingZeros,signalDuration,trailingZeros,measurementDuration,hWater,xSource,ySource,zSource,xRec,yRec,zRec                             