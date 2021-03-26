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















def TankMode(perm=10,fmin=0,fmax=1000,Lx=1.22,Ly=3.66,Lz=0.6,c=1478,walls='rigid',alpha=0,plot=True,pstyle='colored',num=4):
    """
    Parameters
    ----------
    perm:   float, Optional;
            number of permutations to iterate through for each dimension x, y, z. 
            Will end up with arrays for f and mode of size perm**3. Default is 10.
    fmin:   float, Optional;
            Minimum frequency of interest to reduce computation time and limit 
            output array size. Defaults to 0 Hz
    fmax:   float, Optional;
            Maximum frequency of interest to reduce computation time and limit
            output array size. Defaults to 1000 Hz
    Lx:     float, Optional;
            Width of water tank. Defalts as 1.22 m for the BYU Underwater tank. 
            This could be altered when anechoic panels are placed in the tank. 
    Ly:     float, Optional;
            Length of water tank. Defalts as 3.66 m for the BYU Underwater tank. 
            This could be altered when anechoic panels are placed in the tank. 
    Lz:     float, Optional;
            Depth of water in the tank. Defalts as 0.6 m for the BYU Underwater 
            tank. This SHOULD be altered dependent on the current water level
            in the tank. 
    c:      float, Optional;
            Speed of sound in water. This defaults to 1478 m/s following Garrett's
            formula for speed of sound due to Depth, Salinity, and Temperature. 
            This Default is set to the average room temperature of the water and
            assuming near zero salinity over any depth the tank can handle. 
    walls:  string, Optional;
            Choice of hardness of tank walls. This chooses which model to compute. 
            Defaults to walls='rigid' assuming perfectly rigid walls and floor 
            with a pressure release surface.
            ##########################UPDATE!!!!##################################
            walls = 'lossy' assumes lossy walls due to absoptivity alpha. 
            walls = 'multi' assumes lossy walls and requires multiple inputs for
            alpha
    alpha:  float; Optional;
            Absorptivity of the walls if known. Otherwise assumes perfectly rigid. 
    plot:   Boolian; Optional;
            Choose whether or not to plot the EigenFunctions of the natural 
            frequencies in the x-y, x-z, and y-z planes for the first "num" of 
            modes. Default is set as True to plot. False will not plot.  
    pstyle: string; Optional;
            Defaults to 'colored' to plot contourf plots. Can also choose 'line'
            to plot contour line plots. The latter is only recommended when solving
            for very low frequencies. 
    num:    float; Optional;
            Number of modes to plot if Plot = True. Default is set to 4 modes for
            each coordinate plane for a total of 12 plots (if is shows anything)
           
    Returns
    -------
    f:      ndarray of float;
            Ordered Natural frequencies of the tank environment assuming Rigid 
            walls and a pressure release surface.
    mode:   ndarray of int;
            Associated mode numbers of the natural frequencies of the tank. 
    psi:    3D ndarray of float;
            Ordered Eigenfunctions of the natural frequencies. 

    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the natural frequencies and room modes as defined by Garrett
    eq. 13.12 altered by kz in 13.14 for pressure release boundary
    (aka solutions to the eigenfrequencies for a rigid walled tank with pressure
    release water-air interface)
    
    last modified 3/18/2021      
    """
    
    import numpy as np
    #if no input is given for alpha and it defaults to 0, then use rigid wall solution
    if alpha == 0:
        walls = 'rigid'
    #Perfectly rigid wall solutio
    if walls == 'rigid':
        print('Solving for natural frequencies assuming perfectly Rigid walls')
        #rigid wall solution for natural frequencies & Pressure release surface(nz component)
        fN = lambda nx,ny,nz: c/(2)*np.sqrt((nx/Lx)**2 + (ny/Ly)**2 + ((2*nz+1)/(2*Lz))**2)
        #create empty lists to populate
        mode = []
        f = []
        #iterate through the permutations selected for each mode possibility nx,ny,nz
        for nx in range(0,perm):
            for ny in range(0,perm):
                for nz in range(0,perm):
                    #for only values within the chosen bandwidth fmin<= f <=fmax
                    temp = fN(nx,ny,nz)
                    if temp >=fmin:
                        if temp <= fmax:
                            f.append(fN(nx,ny,nz))
                            mode.append([nx,ny,nz])
                    #while fmin <= fN(nx,ny,nz) <= fmax: 
                        #f.append(fN(nx,ny,nz))
                        #mode.append([nx,ny,nz])
        f = np.array(f)
        mode = np.array(mode)
        idxs = np.argsort(f) #order all the frequencies & associated modes in numerical order of freq.
        f = f[idxs]
        mode = mode[idxs]
        print(f'{len(f)} frequencies recorded in range {fmin}<=f<={fmax}')
    
    if walls == 'lossy':
        print('Solving for natural frequencies assuming perfectly lossy walls')
    
    if walls == 'multi':
        print('Solving for natural frequencies assuming perfectly multiple lossy walls')
        
    """
    #Eigen-Function for rectangular tank assuming rigid walls and pressure release water-air interface
    Psi = lambda nx,ny,nz :np.cos(nx*np.pi*x/Lx) * np.cos(ny*np.pi*y/Ly) * np.cos((2*nz-1)*np.pi*z/(2*Lz))
    x = np.linspace(0,Lx)
    y = np.linspace(0,Ly)
    z = np.linspace(0,Lz)
    psi = []
    print('')
    print('calculating EigenFunctions')
    for i in range(len(mode)):
            psi.append(Psi(mode[i,0],mode[i,1],mode[i,2]))
    psi = np.array(psi)
    
    if plot == True:
        print(f'plotting first {num} EigenFunctions')
        import matplotlib.pyplot as plt
        #number of modes we are interested in contour plotting
        start = 0   #zeroeth mode 0, 0, 0 is weird?. 
        #modes and frequencies of interest
        modeint = mode[start:num+start]
        fint = f[start:num+start] 
        
        #plot over x-y plane w/ z = 0
        #create spatial arrays for the 3-dimensions of the tank
        x = np.linspace(0,Lx)
        y = np.linspace(0,Ly)
        z = 0
        x,y = np.meshgrid(x,y)
        for i in range(len(modeint)):
            psi1 = Psi(modeint[i,0],modeint[i,1],modeint[i,2])
            #check to see if mode actually present in this plane, if not, do not plot
            check =np.ones((len(x),len(y)))
            if np.any(psi1 != check) == True:
                fig,ax=plt.subplots(1,1)
                if pstyle == 'line':
                    cb = ax.contour(x,y,psi1,colors='black',linestyles='dashed')
                    ax.clabel(cb,inline=True,fontsize=15)
                else: 
                    cb = ax.contourf(x,y,psi1)
                    fig.colorbar(cb)
                ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)} Hz where Z={z}m')
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Y (m)')
                plt.show()
            else: 
                print('undesired mode not plotted in x-y')
        
        #plot over x-z plane w/ y = 0
        #create spatial arrays for the 3-dimensions of the tank
        x = np.linspace(0,Lx)
        z = np.linspace(0,Lz)
        y = 0
        x,z = np.meshgrid(x,z)
        for i in range(len(modeint)):
            #iterate through calculating each eigenfunction
            psi2 = Psi(modeint[i,0],modeint[i,1],modeint[i,2])
            #check to see if mode actually present in this plane, if not, do not plot
            check =np.ones((len(x),len(z)))
            if np.any(psi2 != check) == True:
                fig,ax=plt.subplots(1,1)
                if pstyle == 'line':
                    cb = ax.contour(x,z,psi2,colors='black',linestyles='dashed')
                    ax.clabel(cb,inline=True,fontsize=15)
                else:
                    cb = ax.contourf(x,z,psi2)
                    fig.colorbar(cb)
                ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)} Hz where Y={y}m')
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Z (m)')
                plt.show()
            else: 
                print('undesired mode not plotted in x-z')        
        
        #plot over y-z plane w/ x = 0
        #create spatial arrays for the 3-dimensions of the tank
        y = np.linspace(0,Ly)
        z = np.linspace(0,Lz)
        x = 0
        y,z = np.meshgrid(y,z)
        for i in range(len(modeint)):
            #iterate through calculating each eigenfunction
            psi3 = Psi(modeint[i,0],modeint[i,1],modeint[i,2])
            #check to see if mode actually present in this plane, if not, do not plot
            check =np.ones((len(y),len(z)))
            if np.any(psi3 != check) == True:
                fig,ax=plt.subplots(1,1)
                if pstyle == 'line':
                    cb = ax.contour(y,z,psi3,colors='black',linestyles='dashed')
                    ax.clabel(cb,inline=True,fontsize=15)
                else:    
                    cb = ax.contourf(y,z,psi3)
                    fig.colorbar(cb)
                ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)} Hz where X={x}m')
                ax.set_xlabel('Y (m)')
                ax.set_ylabel('Z (m)')
                plt.show()
            else: 
                print('undesired mode not plotted in y-z')
    """
    return f, mode#, Psi 
    
    
    
    
    


########################
# How does TL affect the T60 to adjust for underwater or if we do a T10
# Water attenuation?
########################