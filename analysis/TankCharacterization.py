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
    Section 3.4 equations 3.4.1, 3.4.2, & 3.4.3. Where GembaEQ3.4.1 is the 
    Eyring equation that can be found in 661 notes eq 4-2.4.116 calculated from 
    the overall estimated Spatially Averaged Absorption Coefficient. 
    
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
    Aw = alpha_air*A_waterair/S #water absorption coefficient spatially averaged
    Ai = alpha_acrylic*(A_acrylic)/S #acrylic absorption coefficient spatially averaged
    Absorb = Aw + Ai #spatially averaged Absorption Coefficient
    
    #Eyring equation found in 661 notes eq. 4-24.116
    T60 = (24*np.log(10)/c) * (V/(-S*np.log(1-Absorb)))
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





def alpha_prop(f,T=16,S=5,pH=7.7,depth=0.0006):
    """
    Absorption Coefficient from Propagation losses through Sea Water
    
    Parameters
    ----------
    f:      ndarray of float;
            frequency array for bandwidth of interest for freq. dependent absorption
    T:      float, Optional;
            Temperature of the water in Celcius. Defaults to 16 degrees C.
            Effective for -6<T<35 degrees C. 
    S:      float, Optional;
            Salinity of the water in ppt. Effective for 5<S<50 ppt. Defaults to 
            S = 5 ppt
    pH:     float, Optional;
            pH level of the water. Defaults to 7.7 (though this is high relative
            to the test strips and a normal pool). Effective for 7.7<pH<8.3
    depth:  float, Optional;
            depth of water in km. Defaults to 0.0006 km or 0.6m. Which will 
            make that term in the function negligible as basically zero. 
            Effective for 0<z<7 km. 
    
    Returns
    -------
    alpha_prop:     ndarray of float;
                    absorption coefficient if sound (alpha) for propagation 
                    losses through the water. 
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Primarily due to viscous effects above 100kHz (high), but also due to 
    Chemical relaxation of Boric Acid up to a few kHz (low), Chemical 
    relaxation of Magnesium Sulfate up to a few 100kHz (mid). This formulation 
    comes from Ainslie & McColm 1997 - "A simplified formula for viscous and
    chemical absorption in sea water" Published in JASA 103 equation 2 & 3.
    
    Can apply this in propagation models similar to 

    last modified 4/27/2021      
    """
    import numpy as np
    #relaxation frequency for boron
    f1 = 0.78*(S/35)**0.5*np.exp(T/26)
    #relaxation frequency for magnesium
    f2 = 42*np.exp(T/17)
    
    term1 = 0.106*(f1*f**2)/(f**2+f1**2)*np.exp((pH-8)/0.56)
    term2 = 0.52*(1+T/43)*(S/35)*(f2*f**2)/(f**2+f2**2)*np.exp(-depth/6)
    term3 = 0.00049*f**2*np.exp(-(T/27 + depth/17))
    alpha_prop = term1 + term2 + term3 
    return alpha_prop




def T60meas(ht,fs, d = 0.6,c = 1478,rt = 'T60',plot = False,acc=False,alpha_p=0):
    """
    Calculate the T20, T30, or T60 from Backward Schroeder Integration on the 
    measured impulse response hsys.
    
    Parameters
    ----------
    ht:         ndarray of float;
                Measured Impulse Response of the environment. 
    fs:         float;
                Sampling frequency of the impulse reponse. 
    d:          float, Optional;
                depth of water. Defaults to a common 0.6m of water in the tank
    c:          float, Optional;
                speed of sound in water. Defaults to 1478 rounded to nearest whole
                value for any depth of water the tank can have, using Garrett's eq.
                for the speed of sound in water relative to temperature, depth, and 
                salinity for a temparature of 19 degrees C (rough avg. in tank).
    rt:         String, Optional;
                Choose desired Reverb Time (rt) as T10, T20, T30, or T60. Defaults
                to T60. Choosing less than T60 estimates the T60 by assuming linear
                relationship between chosen rt and T60. 
    plot:       boolian, Optional;
                Defaults to False so as to not Plot the 10log(h(t)**2) and the 
                associated Decay Curve. True would plot the two. 
    acc:        boolian, Optional;
                Account for the assumption that the water-air boundary is perfectly
                reflective. Defaults to False to not make this assumption and give
                the overall spatially averaged absorption coefficient. If True, 
                then the spatially averaged absorption coefficient that is returned
                only accounts for the walls and the floor of the enclosure. 
    alpha_p:    float or ndarray of float, Optional;
                Absorption coefficient due to thermoviscous molecular propagation
                losses. Defaults as 0 such that there is no propagation absorption.
                Can use alpha_prop(f,T,S,pH,depth) code to feed in an array of 
                frequency dependent absorption coefficients due to propagation 
                losses through the water. 
                

    Returns
    -------
    T60:            float;
                    Calculated reverberation time (T60) in the tank in seconds.
                    This is calculated using the Through The System (TTS) 
                    response to evaluate reverberation only in the tank. 
                    (i.e. time it takes for the signal in the tank to drop by 
                    60dB)
    alpha_S:        float;
                    Estimated spatially averaged absorption coefficient for the
                    room based on the measured T60 and the Eyring Equation.  
    

    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the measured T60 in the tank. 
    
    last modified 4/27/2021      
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pylab
    from scipy import stats
    params = {'legend.fontsize': 15,
              'figure.figsize': (15, 10),
             'axes.labelsize': 24,
             'axes.titlesize':28,
             'axes.titleweight':'bold',
             'xtick.labelsize':'xx-large',
             'ytick.labelsize':'xx-large',
             'lines.linewidth':2}
    pylab.rcParams.update(params)
    
    #Backward Schroeder Integration
    #Some guidance for this part found here: 
    #https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/room.py
    #the above link provides an alternate method to more generalize this solution 
    T = 1/fs
    ht1 = ht[:int(0.025*len(ht))]
    #ht1 = ht[:int(0.1*len(ht))]
    schroeder = np.empty(len(ht1))
    for n in range(len(ht1)):
        #backwards integration because integrating the rest of the array from a point n
        schroeder[n] = np.sum(ht1[n:]**2) * T
    
    schroeder_dB = 10*np.log10(schroeder) 
    
    if rt == 'T10':
        #determine T30 between -5dB and -15dB of the max value of the decay curve
        init = -5.0
        end = -15.0
        factor = 6.0 #amount to mult. T10 by to extrapolate T60
    if rt == 'T20':
        #determine T20 between -5dB and -25dB of the max value of the decay curve
        init = -5.0
        end = -25.0
        factor = 3.0 #amount to mult. T20 by to extrapolate T60
    if rt == 'T30':
        #determine T30 between -5dB and -35dB of the max value of the decay curve
        init = -5.0
        end = -35.0
        factor = 2.0 #amount to mult. T30 by to extrapolate T60
    if rt == 'T60':
        #determine T30 between -5dB and -35dB of the max value of the decay curve
        init = -5.0
        end = -65.0
        factor = 1.0 #amount to mult. T30 by to extrapolate T60
    
    
    # Linear regression
    #determine the value on the decay curve where it is nearest the init and end
    #values below the maximum of the the decay curve
    sch_init = schroeder_dB[np.abs(schroeder_dB - init).argmin()]
    sch_end = schroeder_dB[np.abs(schroeder_dB - end).argmin()]
    #indices of where the decay curve matches the init and end condition
    init_sample = np.where(schroeder_dB == sch_init)[0][0]
    end_sample = np.where(schroeder_dB == sch_end)[0][0]
    x = np.arange(init_sample, end_sample + 1) / fs
    y = schroeder_dB[init_sample:end_sample + 1]
    slope, intercept = stats.linregress(x, y)[0:2]
    
    """ dont think this works at all. does not get close to measured
    #Reverberation time (T30)
    #convert samples to time and determine the difference
    t_init = init_sample / fs
    t_end = end_sample / fs
    T1 = t_end - t_init
    T1 = factor*T1"""
    
    RT = np.abs(slope/fs)
    T60 = RT *factor
    print('')
    #print(f'{rt} from decay curve slope =',T1,'s')
    print('T60 from decay curve slope =',np.around(T60,decimals=8),'s')
    print('')
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    V = Lx*Ly*d #volume relative to current water depth
    
    if acc == True:
        S = 2*(Ly*d+d*Lx)+Lx*Ly #total enclosed surface area minus air-water surface
    else: 
        S = 2*(Lx*Ly+Ly*d+d*Lx) #total enclosed surface area including air-water
    #solving Eyring equation w/propagation loss found in 661 notes eq 4-2.4.124
    #with default of alpha_p=0 this simplifies to eq. 4-2.4.116
    alpha_S = 1-np.exp(V/S*(8*alpha_p - (24*np.log(10))/(c*T60)))
    
    if plot == True:
        t = np.linspace(0,len(ht1)/fs,len(ht1))*1000 #converted to ms from s
        Level = 10*np.log10(ht1**2)
        plt.figure()
        #plot the IR**2 in dB
        plt.plot(t,Level)
        plt.xlabel('Time (ms)')
        plt.ylabel('Level (dB)')
        plt.grid()
        #plot Decay Curve
        plt.plot(t,schroeder_dB)
        plt.legend([r'$10log[h^{2}(t)]$','Decay Curve'])
        est,_,_ = T60est(d,c)
        plt.title(f'T60meas={np.around(T60*1000,decimals=2)}ms, T60est={np.around(est*1000,decimals=2)}ms & '+r'$\langle \alpha \rangle_{S} =$'+f'{np.around(alpha_S,decimals=4)}')    
    
    return T60, alpha_S





def TankMode(perm=10,fmin=0,fmax=1000,Lx=1.22,Ly=3.66,Lz=0.6,c=1478):
    """
    Determine the rigid wall condition Eigenmodes, Eigenfrequencies, and 
    Eigenfunctions for any frequency range and a chosen number of permutations
    
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
    
           
    Returns
    -------
    f:      ndarray of float;
            Ordered Natural frequencies of the tank environment assuming Rigid 
            walls and a pressure release surface.
    mode:   ndarray of int;
            Associated mode numbers of the natural frequencies of the tank. 

    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the natural frequencies and room modes as defined by Garrett
    eq. 13.12 altered by kz in 13.14 for pressure release boundary
    (aka solutions to the eigenfrequencies for a rigid walled tank with pressure
    release water-air interface)
    
    last modified 4/7/2021      
    """
    
    import numpy as np
    #Perfectly rigid wall solution
    print('Solving for natural frequencies assuming perfectly Rigid walls')
    #rigid wall solution for natural frequencies & Pressure release surface(nz component)
    fN = lambda nx,ny,nz: c/(2)*np.sqrt((nx/Lx)**2 + (ny/Ly)**2 + ((2*nz-1)/(2*Lz))**2)
    #create empty lists to populate
    mode = []
    f = []
    #iterate through the permutations selected for each mode possibility nx,ny,nz
    #nx,ny: 0,1,2,3...
    for nx in range(0,perm):
        for ny in range(0,perm):
            #nz: 1,2,3...
            #Garrett pg.721 "The nz = 0 solution does not exist since constant
            #pressure in the z-direction is not an option that satisfies the 
            #boundary conditions at z=Lz and z=0 simultaneously."
            for nz in range(1,perm):
                #for only values within the chosen bandwidth fmin<= f <=fmax
                temp = fN(nx,ny,nz)
                if temp >=fmin:
                    if temp <= fmax:
                        f.append(fN(nx,ny,nz))
                        mode.append([nx,ny,nz])
                
    f = np.array(f)
    mode = np.array(mode)
    idxs = np.argsort(f) #order all the frequencies & associated modes in numerical order of freq.
    f = f[idxs]
    mode = mode[idxs]
    print(f'{len(f)} frequencies recorded in range {fmin}<=f<={fmax}')
    return f, mode
    
    
    
def TankFunc(x,y,z,f,mode,Lx=1.22,Ly=3.66,Lz=0.6,plot=True,pstyle='colored',num=4):
    """
    Parameters
    ----------
    x:      float or array;
            Single value for a single position in the tank, or array of values
            to iterate over. 
    y:      float or array;
            Single value for a single position in the tank, or array of values
            to iterate over. 
    z:      float or array;
            Single value for a single position in the tank, or array of values
            to iterate over. 
    f:      Ndarray of float;
            Ordered Eigenfrequency array calculated and output by TankMode function
    mode:   Ndarray of float;
            Ordered Eigenmode array calculated and output by TankMode function
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
    plot:   Boolian; Optional;
            Choose whether or not to plot the EigenFunctions of the natural 
            frequencies in the x-y, x-z, and y-z planes for the first "num" of 
            modes. Default is set as True to plot. False will not plot. This only
            plots if len(x) or len(y) or len(z) != 1 
    pstyle: string; Optional;
            Defaults to 'colored' to plot contourf plots. Can also choose 'line'
            to plot contour line plots. The latter is only recommended when solving
            for very low frequencies. 
    num:    float; Optional;
            Number of modes to plot if Plot = True. Default is set to 4 modes for
            each coordinate plane for a total of 12 plots (if is shows anything)
            
    Returns
    -------
    psi:    3D ndarray of float;
            Ordered Eigenfunctions of the natural frequencies. 
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the natural frequencies and room modes as defined by Garrett
    eq. 13.12 altered by kz in 13.14 for pressure release boundary
    (aka solutions to the eigenfrequencies for a rigid walled tank with pressure
    release water-air interface)
    
    last modified 4/7/2021      
    """
    import numpy as np
    #Eigen-Function for rectangular tank assuming rigid walls and pressure release water-air interface
    Psi = lambda nx,ny,nz :np.cos(nx*np.pi*x/Lx) * np.cos(ny*np.pi*y/Ly) * np.cos((2*nz-1)*np.pi*z/(2*Lz))
    
    psi = []
    print('')
    print('calculating EigenFunctions')
    for i in range(len(mode)):
            psi.append(Psi(mode[i,0],mode[i,1],mode[i,2]))
    psi = np.array(psi)

    
    if len(x) > 1:
        ################################################################################
        #The rest of this function is solely for contour plotting the Eigenfunctions psi
        ################################################################################
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
    else:
        #not sure exactly what this plot physically means, but I have it here for now. 
        #might need to change this to default to not plotting if singular value is input.
        if plot == True:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(f,psi)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel(r'$\Psi (r)$')
            plt.title(rf'Eigenfunction $\Psi$({x},{y},{z})')
    
    return psi 



########################
# How does TL affect the T60 to adjust for underwater or if we do a T10
# Water attenuation?
########################