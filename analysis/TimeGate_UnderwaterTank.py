# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:23:48 2020

@author: cvongsaw
"""

###############################################################################
#*****************************************************************************#
#*****************************************************************************#
#*************************** INPUTS VALUES BELOW *****************************#
#*****************************************************************************#
#*****************************************************************************#
###############################################################################

def gatefunc(IR,fs,tgate,tb4=0.1):
    """
    Compute the Sound Speed of the Water based on the Depth, Temperature, and 
    Salinity of the water according to three well known models. Garrett, Medwin 
    & Kuperman or Wilson. 
    
    Parameters
    ----------
    IR:     ndarray;
            Impulse Response or time domain signal. 
    fs:     float;
            Sampling frequency of the input time domain signal. Measured in Hz.
    tgate:  float;
            Amount of time in seconds from the beginning of the input IR signal
            in which the reflection of interest is arriving that needs to be
            timegated out of the signal.
    tb4:    float, optional;
            Defaults to 0.1 ms. This is the time before the reflection that the 
            timegating should start to cut off any buildup to the reflected signal. 
            This should also ideally be after the initial direct signal. 

    Returns
    -------
    IRgate: ndarray;
            Time-gated array of the input signal.  
    Notes
    -----
    Author: Cameron Vongsawad
    
    
    Last Modified: 4/1/2021
    """
    import numpy as np
    Nb4 = tb4/1000 *fs #convert to seconds and then samples before gating
    #where to start the time-gating or cutting off the signal to zero
    fin = int(tgate*fs-Nb4) 
    start = 0 #start the time gating allowing everything from the beginning of ht
    #convert time length to samples to determine the finish cutoff of ht
    #fin = int(tgate*fs*percent)
    #cut off the IR before the first reflection being index "fin"
    IRgate = np.zeros(len(IR))
    IRgate[start:fin] = IR[start:fin] #replace up to gate with original
    tbuff = tb4/2 #buffer value to determine where to apply the hanning window.
    damp = int(tbuff/1000*fs) #0.05ms of damping converted to samples
    #apply hanning window to portion of the array following original cutoff
    #this allows for the signal to more gradually ramp down to zeros. 
    IRgate[fin:fin+damp] = IR[fin:fin+damp]*np.hanning(damp)
    #repopulate first half of that damping data keeping original array information
    IRgate[fin:int(fin+damp/2)] = IR[fin:int(fin+damp/2)]
    
    return IRgate


def uwsoundspeed(D=0.2,T=16.0,S=0.03,model='Garrett'):
    """
    Compute the Sound Speed of the Water based on the Depth, Temperature, and 
    Salinity of the water according to three well known models. Garrett, Medwin 
    & Kuperman or Wilson. 
    
    Parameters
    ----------
    #Water Characteristics in the Tank#
    D:  float, optional;
        water depth (m) where 0<= D <=1000m
    T:  float, optional;
        temperature in Celcius where -2<= T <=24.5 
    S:  float, optional;
        salinity where 0.030<= S <=0.042 grams salt per kg H20 (aka parts per
        thousand = ppt)
    Returns
    -------
    c:  float;
        speed of sound in water for the specified depth, temperature and salinity   
    
    Notes
    -----
    Author: Cameron Vongsawad

    Last Modified: 4/1/2021
    """
    ###############################################################################
    ######## sound speed (m/s), Cite Garrett valid w/in +-0.2m/s ##################
    ###### appears to be accurate w/in 0.000969% of wiki value @20C ###############
    ########## effects of depth is negligible in the tank limits ##################
    ########### EQ 11.26 pg 619 Garrett "Understanding Acoustics" #################
    ###############################################################################
    if model == 'Garrett' or 'garrett' or None:
        c = 1493 + 3*(T-10) - 0.006*(T-10)**2 - 0.04*(T-18)**2 + 1.2*(S-35)- 0.01*(T-18)*(S-35) + D/61
    
    #medwin & Kuperman Encyclo. of Ocean Sciences 2nd ed. 2001
    if model == 'Kuperman' or 'Medwin' or 'MedwinKuperman':
        c = 1449.2 + 4.6*T - 0.055*T**2 + 0.00029*T**3 + (1.34 - 0.010*T)*(S - 35) + 0.016*D

    #Christ, WenliSr., The ROV Manual 2nd ed. 2014 from simplified Wilson's 1960
    #S is in PSU which is basically equivalent to ppt
    if model == 'Wilson':   
        c = 1449 + 4.6*T - 0.055*T**2 + 0.0003*T**3 + 1.39*(S - 35) + 0.017*D
    return c


def gateValue(AEgir_pose, Ran_pose, D, c=1478, Coordinate='tank',Print='True'):
    import numpy as np
    """
    Compute the first bounce reverberations of the BYU Hydroacoustics lab tank
    in order to timegate signals. Assumes a rectangular volume. 
    
    Parameters
    ----------
    AEgir_Pose: tuple; 
                AEgir TCP position (x,y,z)
    Ran_pose:   tuple;
                Ran TCP position (x,y,z,v) if 'robot' frame or (x,y,z) if 'tank' frame
                where v is the vention position (7th axis extender)   
   
    #Water Characteristics in the Tank#
    D:  float, optional;
        water depth (m) where 0<= D <=1000m
    
    Coordinate: string, optional;
                Choose if cordinate system is robot frame or tank frame
                Standard is tank frame "tank"
                or robot frame inputing each robot + vention positioning "robot" 
    Print: string, optional;
                Choose if you want the function to print a bunch of numbers


    Returns
    -------
    tshort: float;
            shortest time for a single reflection in seconds. 
    tside:  float
            shortest time for single reflection of side wall reflections only 
            but still allowing potential for seabed and surface reflections.
    directpath: float;
                distance of direct path from hydrophone to hydrophone
    
    prints values of: 
        AEgir and Ran positions
        direct sound "tdirect"
        Single bounce times:
            bottom "tb"
            H2O-O2 "tt"
            Side 1 "ts1"
            Side 2 "ts2"
            Front wall "tfront"
            Back wall "tback"
            
    Notes
    -----
    Author: Cameron Vongsawad
    
    This code only allows for a single tuple of length len(A)=3 and len(R)=3 or 4
    
    
    Last Modified: 11/17/2020
    
    
    Times printed in "ms" (milliseconds), however programmed values in seconds
    """          
        ###############################################################################
        ################### function to determine time of flight for ##################
        ############## ray paths knowing tank fram positions ##########################
        ###############################################################################
    def pathtime(XA,YA,ZA,XR,YR,ZR):
        """
        XA, YA, ZA  : float, cartesian coordinates of AEgir
        XR, YR, ZR  : float, cartesian coordinates of Ran
            
        """
        ###############################################################################
        ######################### main code for direct path time ######################
        ###############################################################################
        directpath = np.sqrt((XA-XR)**2+(YA-YR)**2+(ZA-ZR)**2)      #direct distance eq
        tdirect = (directpath)/c
        
        ###############################################################################
        ######################### main code for bottom bounce #########################
        ####################### determined through geometries #########################
        ###############################################################################
        range_b = np.sqrt(directpath**2 - np.abs(ZA-ZR)**2)         #r bottom z-y plane
        rzA = ZA/np.sin(np.arctan((ZA+ZR)/range_b))
        rzR = ZR/np.sin(np.arctan((ZA+ZR)/range_b))
        tb = (rzA + rzR)/c                                          #bottom bounce time
        
        ###############################################################################
        ######################### main code for top bounce ############################
        ####################### determined through geometries #########################
        ###############################################################################
        ZAt = D - ZA                           #translate to looking from water surface
        ZRt = D - ZR
        range_t = np.sqrt(directpath**2 - np.abs(ZAt-ZRt)**2)       #r top  z-y plane
        rzAt = ZAt/np.sin(np.arctan((ZAt+ZRt)/range_t))
        rzRt = ZRt/np.sin(np.arctan((ZAt+ZRt)/range_t))
        tt = (rzAt + rzRt)/c                                        #top bounce time
        
        ###############################################################################
        ######################### main code for side1 x=0 bounce ######################
        ####################### determined through geometries #########################
        ###############################################################################
        range_s1 = np.sqrt(directpath**2 - np.abs(XA-XR)**2)        #r x=0 x-y plane
        rxAs1 = XA/np.sin(np.arctan((XA+XR)/range_s1))
        rxRs1 = XR/np.sin(np.arctan((XA+XR)/range_s1))
        ts1 = (rxAs1 + rxRs1)/c                                     #x=0 bounce time
        
        ###############################################################################
        ######################### main code for side2 X=X bounce ######################
        ####################### determined through geometries #########################
        ###############################################################################
        Xmax = 1.22
        XAs2 = Xmax - XA                       #translate to looking from water surface
        XRs2 = Xmax - XR
        range_s2 = np.sqrt(directpath**2 - np.abs(XAs2-XRs2)**2)    #r x=x x-y plane
        rxAs2 = XAs2/np.sin(np.arctan((XAs2+XRs2)/range_s2))
        rxRs2 = XRs2/np.sin(np.arctan((XAs2+XRs2)/range_s2))
        ts2 = (rxAs2 + rxRs2)/c                                     #x=x bounce time
          
        ###############################################################################
        ################## main code for "front" (North) wall bounce 1 y=0 ############
        ################## using the method of images                      ############
        ###############################################################################
        range_front = np.sqrt((XA-XR)**2+(YA-(-YR))**2+(ZA-ZR)**2)  #direct image
        tfront = (range_front)/c                                    #front wall time
        
        ###############################################################################
        ################## main code for "back" (South) wall bounce 2 y=y##############
        ################## using the method of images                      ############
        ###############################################################################
        Ymax = 3.66
        range_back = np.sqrt((XA-XR)**2+(YA-(YR+Ymax))**2+(ZA-ZR)**2)   #direct image
        tback = (range_back)/c                                      #front wall time
        
        
        
        t = (tb,tt,ts1,ts2,tfront,tback)
        tshort = min(t)
        tside = (ts1,ts2,tfront,tback)
        tside = min(tside)
        if Print:
            print('')
            print('AEgir(source) & Ran(Receiver) tank frame coordinates:')
            print('(XA,YA,ZA)=',(XA,YA,ZA))
            print('(XR,YR,ZR)=',(XR,YR,ZR))
            print('')
            print('Single Bounce reverberation times to receiver:')
            print('direct sound t=', tdirect*10**3,'ms')
            print('bottom bounce t=', tb*10**3,'ms')
            print('H20-O2 bounce t=', tt*10**3,'ms')
            print('Side 1 bounce t=', ts1*10**3,'ms')
            print('Side 2 bounce t=', ts2*10**3,'ms')
            print('Front Wall bounce t=', tfront*10**3,'ms')
            print('Back Wall bounce t=', tback*10**3,'ms')
            print('tshort=',tshort,'s')
            print('')
        return tshort,tside,tdirect,directpath     
        
              
        #### for tank frame coordinates, no need to translate coordinates #############
    if Coordinate == 'tank': 
        ###############################################################################
        ## Hydrophone Locations (Insert AEgir & Ran Tank coordinates (X,Y,Z) in m) ####
        ###############################################################################
                        #"AEgir" Tank Frame position (X,Y,Z) TCP TC4038  
        XA   = AEgir_pose[0]    
        YA   = AEgir_pose[1] 
        ZA   = AEgir_pose[2]
                        #"Ran" Tank Frame position (X,Y,Z) TCP TC4034 
        XR   = Ran_pose[0]           
        YR   = Ran_pose[1]            
        ZR   = Ran_pose[2]
        
        tshort,tside,tdirect,directpath = pathtime(XA,YA,ZA,XR,YR,ZR)
    
        
        #### for robot frame coordinates, must translate to tank frame first ##########
    elif Coordinate == "robot":      
        ###############################################################################
        ########## TCP Locations (insert current TCP locations in mm)##################
        #### This is in correlation with default settings for the end connector ####### 
        ###############################################################################
                        #"AEgir" position TCP TC4038  
        TCPxA   = AEgir_pose[0]    
        TCPyA   = AEgir_pose[1] 
        TCPzA   = AEgir_pose[2]
                        #"Ran" position TCP TC4034 
        TCPxR   = Ran_pose[0]           
        TCPyR   = Ran_pose[1]            
        TCPzR   = Ran_pose[2]
        TCPvR   = Ran_pose[3]    #Vention 7th axis positioning adjustment for y direction     
        
        ###############################################################################
        ######## Tank Frame Locations (insert current tank locations in mm)############
        ######## these are directly measured values. must comment out future ##########
        ######## translation of positioning if used. OR translation trumps this #######
        ###############################################################################
                        
        #convert mm positioning to m
        TCPxA = TCPxA/1000 
        TCPyA = TCPyA/1000
        TCPzA = TCPzA/1000
        TCPxR = TCPxR/1000      
        TCPyR = TCPyR/1000           
        TCPzR = TCPzR/1000
        TCPvR = TCPvR/1000

        ###############################################################################
        ## translating TCP position to tank coordinate positions (/1000 for mm => m) ##
        ## Home position used for conversion w/end connector settings of both #########
        ## "AEgir" and "Ran" measured in mm initially and then later converted to m ###
        # Home position in the tank frame for "Ran" should be measured at VR = 0 ######
        ###############################################################################
        XA_TCP_home = 392.69/1000
        YA_TCP_home = 288.83/1000
        ZA_TCP_home = -59.54/1000
        XA_tank_home = 100/1000  
        YA_tank_home = 2984/1000 
        ZA_tank_home = 901/1000  
        
        XR_TCP_home = 1291.6/1000 
        YR_TCP_home = 132.98/1000
        ZR_TCP_home = -190.91/1000
        VR_TCP_home = 1404.7
        XR_tank_home = 1010/1000
        YR_tank_home = 541/1000
        ZR_tank_home = 721/1000
        VR_tank_home = YR_tank_home
        
        XA =  (XA_tank_home + (-XA_TCP_home + TCPxA) )
        YA =  (YA_tank_home + (-YA_TCP_home + TCPyA) )
        ZA =  (ZA_tank_home + (-ZA_TCP_home + TCPzA) )
        XR =  (XR_tank_home + (-XR_TCP_home + TCPxR) )
        YR =  (YR_tank_home + (-YR_TCP_home + TCPyR) -TCPvR ) #adjusted for Vention pos
        ZR =  (ZR_tank_home + (-ZR_TCP_home + TCPzR) )
        
        tshort,tside = pathtime(XA,YA,ZA,XR,YR,ZR)

    return tshort,tside,tdirect,directpath
