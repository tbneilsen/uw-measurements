# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:23:48 2020

@author: cvong
"""
import numpy as np
###############################################################################
#*****************************************************************************#
#*****************************************************************************#
#*************************** INPUTS VALUES BELOW *****************************#
#*****************************************************************************#
#*****************************************************************************#
###############################################################################


###############################################################################
########## water characteristics ##############################################
###############################################################################
D = 0.618               #water depth (m) where 0<= D <=1000m
T = 16.0                #temperature in celcius where -2<= T <=24.5 
S = 0.03                #salinity where 0.030<= S <=0.042 grams salt per kg H20

#receive = TC4034 on stationary rod taped to pvc cross
#source = TC4038 on stationary UR10e              
###############################################################################
########## TCP Locations (insert current TCP locations in mm)##################
#### This is in correlation with settings for the end connector found in ######
###### Test 1 Precision and Repeatability Program on the teach pendant ########  
###############################################################################
                #Source TCP TC4038  
TCPxs = -748.2    
TCPys = -799.2
TCPzs = -1141.8
                #Receiver TCP TC4034 
TCPxr = -922           
TCPyr = 500            
TCPzr = -198       

###############################################################################
######## Tank Frame Locations (insert current tank locations in mm)############
######## these are directly measured values. must comment out future ##########
######## translation of positioning if used. OR translation trumps this #######
###############################################################################
                 #UR10e 1 w/ TC4038 Source in tank frame
XS = 540.0              #initial x pos. to translate from robot to tank (mm)
YS = 1967.0             #initial y pos. to translate from robot to tank (mm)
ZS = 330.0              #initial z pos. to translate from robot to tank (mm)
                 #UR10e 2 w/ TC4034 Receiver in tank frame
XR = 510.0              #initial x pos. to translate from robot to tank (mm)
YR = 1910.0             #initial y pos. to translate from robot to tank (mm)
ZR = 195.0              #initial z pos. to translate from robot to tank (mm)
                 #convert mm to m
XS = XS/1000
YS = YS/1000
ZS = ZS/1000    
XR = XR/1000
YR = YR/1000
ZR = ZR/1000
TCPxs = TCPxs/1000 
TCPys = TCPys/1000
TCPzs = TCPzs/1000
TCPxr = TCPxr/1000      
TCPyr = TCPyr/1000           
TCPzr = TCPzr/1000
###############################################################################
#*****************************************************************************#
#*****************************************************************************#
#************************** STOP! NO MORE INPUTS *****************************#
#*****************************************************************************#
#*****************************************************************************#
###############################################################################


###############################################################################
######## sound speed (m/s), Cite Garrett valid w/in +-0.2m/s ##################
###### appears to be accurate w/in 0.000969% of wiki value @20C ###############
###############################################################################
c = 1493 + 3*(T-10) - 0.006*(T-10)**2 - 0.04*(T-18)**2 + 1.2*(S-35)
- 0.01*(T-18)*(S-35) + D/61

###############################################################################
## translating TCP position to tank coordinate positions (/1000 for mm => m) ##
## Home position used for conversion w/end connector settings of both #########
## source and receiver measured in mm initially and then later converted to m #
###############################################################################
XS_TCP_home = 208.8/1000
YS_TCP_home = -477.6/1000
ZS_TCP_home = -565.9/1000
XS_tank_home = 147/1000
YS_tank_home = 2955/1000
ZS_tank_home = 790/1000
"""
XR_TCP_home = 
YR_TCP_home = 
ZR_TCP_home = 
XR_tank_home = 
YR_tank_home = 
ZR_tank_home = 
"""
XS =  (XS_tank_home + ( YS_TCP_home - TCPys) )
YS =  (YS_tank_home + (-XS_TCP_home + TCPxs) )
ZS =  (ZS_tank_home + (-ZS_TCP_home + TCPzs) )
"""
XR =  (XR_tank_home + (YR_TCP_home - TCPyr))
YR =  (YR_tank_home + (-XR_TCP_home + TCPxr))
ZR =  (ZR_tank_home + (-ZR_TCP_home + TCPzr))
"""
print('source and receiver tank frame coordinates:')
print('(XS,XY,XZ)=',(XS,YS,ZS))
print('(XR,XR,XR)=',(XR,YR,ZR))
print('')
print('Single Bounce reverberation times to receiver:')

###############################################################################
######################### main code for bottom bounce #########################
###############################################################################
directpath = np.sqrt((XS-XR)**2+(YS-YR)**2+(ZS-ZR)**2)      #direct distance eq
range_b = np.sqrt(directpath**2 - np.abs(ZS-ZR)**2)         #r bottom z-y plane
rzS = ZS/np.sin(np.arctan((ZS+ZR)/range_b))
rzR = ZR/np.sin(np.arctan((ZS+ZR)/range_b))
tb = (rzS + rzR)/c                                          #bottom bounce time
print('bottom bounce t=', tb*10**3,'millisec')

###############################################################################
######################### main code for bottom bounce #########################
###############################################################################
ZSt = D - ZS                           #translate to looking from water surface
ZRt = D - ZR
range_t = np.sqrt(directpath**2 - np.abs(ZSt-ZRt)**2)       #r top  z-y plane
rzSt = ZSt/np.sin(np.arctan((ZSt+ZRt)/range_t))
rzRt = ZRt/np.sin(np.arctan((ZSt+ZRt)/range_t))
tt = (rzSt + rzRt)/c                                        #top bounce time
print('H20-O2 bounce t=', tt*10**3,'millisec')

###############################################################################
######################### main code for side1 x=0 bounce ######################
###############################################################################
range_s1 = np.sqrt(directpath**2 - np.abs(XS-XR)**2)        #r x=0 x-y plane
rxSs1 = XS/np.sin(np.arctan((XS+XR)/range_s1))
rxRs1 = XR/np.sin(np.arctan((XS+XR)/range_s1))
ts1 = (rxSs1 + rxRs1)/c                                     #x=0 bounce time
print('Side 1 bounce t=', ts1*10**3,'millisec')

###############################################################################
######################### main code for side2 X=X bounce ######################
###############################################################################
Xmax = 1.22
XSs2 = Xmax - XS                       #translate to looking from water surface
XRs2 = Xmax - XR
range_s2 = np.sqrt(directpath**2 - np.abs(XSs2-XRs2)**2)    #r x=x x-y plane
rxSs2 = XSs2/np.sin(np.arctan((XSs2+XRs2)/range_s2))
rxRs2 = XRs2/np.sin(np.arctan((XSs2+XRs2)/range_s2))
ts2 = (rxSs2 + rxRs2)/c                                     #x=x bounce time
print('Side 2 bounce t=', ts2*10**3,'millisec')
  
###############################################################################
################## main code for front wall bounce 1 y=0 ######################
###############################################################################
tfront = 0
print('Front Wall bounce t=', tfront*10**3,'millisec')

###############################################################################
################## main code for back wall bounce 2 y=y #######################
###############################################################################
Ymax = 3.66
tback = 0
print('Back Wall bounce t=', tback*10**3,'millisec')
