# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:50:24 2020

Code developed for basic characterization of the tank environment. Including
T60 estimation, schroeder frequency estimation, trailing zeros

@author: cvongsaw
"""


def T60est(d,c = 1478,zw = 1.5e6,zi= 3.26e6):
    """
    Parameters
    ----------
    d:      float;
            depth of water
    c:      float, Optional;
            speed of sound in water. Defaults to 1478m/s (rounded to nearest 
            whole value for any depth of water the tank can have), using 
            Garrett's eq. for the speed of sound in water relative to 
            temperature, depth, and salinity for a temparature of 19 degrees C 
            (rough avg. in tank).
    zw:     float, Optional;
            Acoustic impedance of water is generally accepted as 1.5E6 Ns/m**-3
            and defaults to this value.
    zi:     float, Optional;
            Acoustic impdeance of side walls. Defaults to acoustic impdedance 
            of acrylic, accepted as 3.26E6 Ns/m**-3 from the following source:
                https://www.ndt.net/links/proper.htm 
                
    Returns
    -------
    T60:            float;
                    Estimate of the reverberation time (T60) in seconds. 
                    i.e. time it takes for the signal to drop by 60dB
    sigL:           float;
                    minimum excitation signal length (s) required by T60 based on Gemba
                    recommendation for 5-10x length of T60. This gives 10x. 
    fschroeder:     float;
                    Schroeder Frequency (Hz). The lowest frequency of interest in which 
                    the tank is large.
    

    Notes
    -----
    Author: Cameron Vongsawad
    
    Comes from Gemba 2014 dissertation ("Characterization of underwater acoustic 
    sources recorded in reverberant environments with application to scuba...")
    Section 3.4 equations 3.4.1, 3.4.2, & 3.4.3
    
    This can then be used to determine the length of the excitation signal (which
    must be 5-10x longer than the T60) 
    
    This is all relative to the depth of the water, and if the boundary 
    impedance is altered from the standard tank. 
    
    Changed order of output
    
    last modified 2/22/2021      
    """
    import numpy as np
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    V = Lx*Ly*d #volume relative to current water depth
    A_acrylic = Lx*Ly + 2*Ly*d + 2*Lx*d #total surface area of acrylic boundaries
    A_waterair = Lx*Ly #total surface area of water-air boundary
    S = A_acrylic +A_waterair #total surface area of (semi)absorptive boundaries
    
    #estimate absorption coefficients for boundaries
    za = 415 #accepted impedance of air in Ns/m**-3
    alpha_acrylic = 1-np.abs((zw-zi)/(zw+zi))
    alpha_air = 1-np.abs((zw-za)/(zw+za))
    
    #Sum of alpha*A/S found in eq. 3.4.1 of Gemba
    #Absorption can be more thoroughly estimated using Physcs 661 notes. 
    Aw = alpha_air*A_waterair/S #water absorption
    Ai = alpha_acrylic*(A_acrylic)/S #acrylic absorption
    Absorb = Aw + Ai
    
    T60 = (24*np.log(10)/c) * (V/(-S*np.log10(1-Absorb)))
    fschroeder = 0.6*np.sqrt(c**3*T60/V)
    signal_length = 10*T60
    sigL = signal_length
    
    #if desired to compare with a simpler room estimation found in 461 notes?
    #T60ng = np.log(10**6)*4*V/(c*Absorb)
    
    return T60, sigL, fschroeder


def trailzeros(RT,sigL,fstart,fstop,f = None,R = 60,sig = 'lin') :
    """
    *****Super not sure if this is working because linear and exponential dont 
    give differing results. and tstop does not seem like it is calculated 
    correctly since when f = None should cause it to give l as solution********
    
    
    Parameters
    ----------
    RT:     float;
            Reverberation time (s). Defaults to the T60, but can be altered by 
            changing the following parameter R. T60 estimate can be determined
            by the depth of the water using the function T60est
    sigL:   float;
            Generated Signal length (s). The min. length can be found in 
            T60est(d) and should be 5-10x that of the estimated T60
    fstart: float;
            Start frequency (Hz) of the chirp signal
    fstop:  float;
            Stop frequency (Hz) of the chirp signal
    f:      float;
            Target Frequency of interest within the chirped signal, often the 
            highest frequency and therefore defaults as None which makes
            f = fstop. Chosen as the highest frequency of interest.
    R:      float, Optional;
            Defaults to 60dB as the dynamic range for the reverberation time 
            (T60), but can be change to a T15, T25, etc. 
    sig:    string, Optional;
            Signal type. Either 'lin' for linear or 'exp' exponential chirp. 
            Defaults to 'lin' chirped signal. 
    
    Returns
    -------
    tstop:  float;
            Trailing zeros necessary (stop margin, or stop gap)
    tls:    float;
            Total length of signal and trailing zeros recommended. 

    Notes
    -----
    Author: Cameron Vongsawad
    
    Trailing Zeros estimate from Muller-Trapet JASA 2020 paper based on RT. 
    
    Changed order of input to align better with T60est()
    
    Last Modified: 2/22/2021
            
"""
    #tf = time in the sweep, when certain frequency f is played
    #D = dynamic range for RT found by D = 20dB + R where R is the reverb time 
    # level decrease (where for a T60 will be R = 60 and D = 80) 
    import numpy as np 
    
    if f == None: 
            f = fstop
    if sig == 'lin' or 'linear':
        #eq 21 time in the sweep when the target frequency occurs (sigL = min. 
        #actual sweep len)
        tlin = sigL*(f-fstart)/(fstop-fstart)
        #eq 20 determining the total signal/recording length duration from 
        #dynamic range and the estimated RT
        l = tlin + (20 + R)/60*RT
        #eq 18 determine the stop margine or time of trailing zeros for the 
        #signal (l = total signal duration, sigL = time of sweep design)
        tstop = l - sigL 
            
    if sig == 'exp' or 'exponential':
        #eq 24 time in the sweep when the target frequency occurs (sigL = min. 
        #actual sweep len) for exponential chirps
        texp = sigL* np.log(f/fstart)/np.log(fstop/fstart)
        #eq 20 determining the total signal/recording length duration from 
        #dynamic range and the estimated RT
        l = texp + (20 + R)/60*RT
        #eq 25 determine the stop margine or time of trailing zeros for the 
        #signal (l = total signal duration, sigL = time of sweep design)
        tstop = ((20+R)/60*texp*np.log(fstop/fstart) - l *np.log(fstop/f))/np.log(f/fstart)
    
    tls = l
    print(sigL)
    print(tstop)
    print(tls)
    return tstop, tls






def T60meas(rec, hsys, d = 0.6,c = 1478,zw = 1.5e6,zi= 3.26e6):
    """
    Parameters
    ----------
    rec:    ndarray of float;
            recorded signal to determine the T60 from. 
    hsys:   ndarray of float;
            Through The System (TTS) response from ESAUResponse SysResponse. 
            This is used to take out the system response and any electrical noise. 
    d:      float;
            depth of water
    c:      float, Optional;
            speed of sound in water. Defaults to 1478 rounded to nearest whole
            value for any depth of water the tank can have, using Garrett's eq.
            for the speed of sound in water relative to temperature, depth, and 
            salinity for a temparature of 19 degrees C (rough avg. in tank).
    zw:     float, Optional;
            Acoustic impedance of water is generally accepted as 1.5E6 Ns/m**-3
            and defaults to this value.
    zi:     float, Optional;
            Acoustic impdeance of side walls. Defaults to acoustic impdedance 
            of acrylic, accepted as 3.26E6 Ns/m**-3 from the following source:
                https://www.ndt.net/links/proper.htm 
                
    Returns
    -------
    T60:            float;
                    Calculated reverberation time (T60) in the tank in seconds.
                    This is calculated using the Through The System (TTS) 
                    response to evaluate reverberation only in the tank. 
                    (i.e. time it takes for the signal in the tank to drop by 
                    60dB)
    

    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the measured T60 in the tank. 
    
    last modified 2/22/2021      
    """
    import numpy as np
    
    
    
    
    return T60



















"""
room modes as defined by Garrett
"""
import numpy as np

D = 0.92                #water depth (m) where 0<= D <=1000m
T = 19.0                #temperature in celcius where -2<= T <=24.5 
S = 0.03                 #salinity where 0.030<= S <=0.042 grams salt per kg water              

#sound speed (m/s), Cite Garrett valid w/in +-0.2m/s
c = 1493 + 3*(T-10) - 0.006*(T-10)**2 - 0.04*(T-18)**2 + 1.2*(S-35) - 0.01*(T-18)*(S-35) + D/61

print('c =', c, 'm/s')
Lx = 1.22   #m
Ly = 3.66   #m
Lz = D      #depth in m


n = 10
mode = []
f = []
#eq. 13.12 Garrett altered by kz in 13.14 for pressure release boundary
def TankMode(nx,ny,nz,Lx=Lx,Ly=Ly,Lz=Lz,c=c):
    return c/(2)*((nx/Lx)**2 + (ny/Ly)**2 + ((2*nz-1)/(2*Lz))**2)**(1/2)

for nx in range(0,n):
    for ny in range(0,n):
        for nz in range(0,n):
            if TankMode(nx,ny,nz) <= 1000: 
                f.append(TankMode(nx,ny,nz))
                mode.append([nx,ny,nz])
            while 100000 <= TankMode(nx,ny,nz) <= 300000: 
                f.append(TankMode(nx,ny,nz))
                mode.append([nx,ny,nz])

f = np.array(f)
mode = np.array(mode)
idxs = np.argsort(f)
f = f[idxs]
mode = mode[idxs]

print('# freq recorded in range',len(f))
print('resonant freqs(Hz) =', f[1:6])#,f[500000:500010])
print(mode[1:6])#,mode[500000:500010])



########################
# How does TL affect the T60 to adjust for underwater or if we do a T10
# Water attenuation?
########################






