# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 19:50:24 2020

Code developed for basic characterization of the tank environment. 
Including:
Fractional octave band filtering
T60 estimation and schroeder frequency estimation
Minimum trailing zeros estimation
Propagation absorption coefficient
Applet for determining time bounds for T60 measurement
T60 measured from h(t)
Absorption coefficient of walls in tank from T60 
Absorption coefficient of material added to the tank
Tank eigenmodes
Tank eigenfunctions
A finite-impedance boundary modal model following Pierce




@author: cvongsaw
"""
def OctaveFilter(data,f0,f1,fs,frac = 1,order = 5,exact = True):
    """
    
    Parameters
    ----------
    data:       Ndarray;
                Sampled data that covers some bandwidth.
    f0:         float;
                Low-end frequency (Hz) of the desired bandwidth
    f1:         float;
                High-end frequency (Hz) of the desired bandwidth
    fs:         float;
                Sampling frequency of the data
    frac:       float, Optional;
                Bandwidth fraction. Examples: 1/3-octave frac=3, 1-octave frac=1 
                (Default), 2/3-octave frac=3/2.
    order:      Int, Optional;
                Order of the filter. Defaults to 5. 
    exact:      boolean;
                Gives option to use IEC standard for octave ratio (10**(3/10))
                or generally accepted standard of 2. Default is True. Set exact
                to False if factor of 2 is desired. 
    
    Returns
    -------
    filt_data:  Ndarray;
                2-d array of the bandpass filtered data. Row dimensions = same 
                dimensions as mid_bands. Each row is the data for a given band. 
                The column dimensions are the filtered data. Ex) filt_data[0,:] 
                would be all of the data for the first mid-band frequency. 
                
    mid_bands:  Ndarray of float;
                Array of octave or fractional octave frequencies
                Note: center frequencies are based on IEC standard 61260-1
                found in equation 1 in section 5.2.1. This code defaults to the 
                octave ratio 10**(3/10) as opposed to the standard ratio of 2. 
    

    Notes
    -----
    Author: Corey Dobbs
    
    Apply a bandpass filter to data in order to obtain an average over an 
    octave or fractional octave band centered at the middle frequencies output
    in mid_bands. 
    
    References:
    https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    
    https://github.com/jmrplens/PyOctaveBand/blob/
    43e65e6cfc50d0b079383fee7ba0693cd645c350/PyOctaveBand.py#L14
    
    TDOTOspec.m by Dr. Kent Gee at BYU, found in BYU Acoustics  
    under General Signal Processing/src/Analyzing Spectra
    https://git.physics.byu.edu/acoustics
    
    Dr. Gee's code included this note:
    BUTTER is based on a bilinear transformation, as suggested in
    ANSI standard.  From oct3dsgn function by Christophe Couvreur, Faculte 
    Polytechnique de Mons (Belgium)

    
    last modified 9/1/2021      
    """
    import numpy as np
    import math
    import scipy.signal as sig
    
    
    #Generate Frequency Array
    if exact == True:    
        G = 10**(3/10) #octave frequency ratio
        #based on IEC standard 61260-1 found in equation 1 in section 5.2.1. 
    elif exact == False:
        G = 2
        #generally accepted octave frequency ratio
    fr = 1000     #reference frequency
    
    
    # Get the initial mid-band frequency
    #According to IEC standard 61260-1 section 5.4
    if frac % 2 == 0: #Even frac
        x_init = math.ceil(frac*np.log(f0/fr)/np.log(G) - 1/2)
        x_final = math.floor(frac*np.log(f1/fr)/np.log(G) - 1/2)
    else: #Odd frac
        x_init = math.ceil(frac*np.log(f0/fr)/np.log(G))
        x_final = math.floor(frac*np.log(f1/fr)/np.log(G))
    
    x = np.arange(x_init,x_final + 1)
    
    
    #Get mid-band frequencies and limits
    if frac % 2 != 0: #Odd frac
        mid_bands = fr*G**(x/frac)
    else: #Even frac
        mid_bands = fr*G**((2*x+1)/(2*frac))
        
    #Get frequency band limits
    #References codes by Kent Gee and Christophe Couvreur
    upper_limits = mid_bands*G**(1/(2*frac)) #low ends of filter
    lower_limits = mid_bands/G**(1/(2*frac)) #high ends of filter   
    Qr = mid_bands/(upper_limits - lower_limits)
    Qd = np.pi/2/frac/np.sin(np.pi/2/frac)*Qr
    alpha = (1 + np.sqrt(1+4*Qd**2))/2/Qd
    
    
    
    #Zero mean
    data = data - np.mean(data)
    
    #Window, and rescaling
    w = np.hanning(len(data))
    data = data*w/np.sqrt(np.mean(w**2))

    
    #Use a butterworth filter on the data according to the fractional octave bands
    for i in range(len(mid_bands)):
        
        #Use a decimation factor to keep the sampling frequency within 
        #reasonable limits.
        
        if mid_bands[i] < fs/20: #factor of 20 suggested as threshold for decimation
            deci_rat = np.ceil(fs/mid_bands[i]/20)  #Decimation factor
            #decdata = sig.decimate(sig.decimate(data,10),2)
        else:
            deci_rat = 1
        
        fsdec = fs/deci_rat #Decimated sampling rate
        
        
        W1 = mid_bands[i]/(fsdec/2)/alpha[i]
        W2 = mid_bands[i]/(fsdec/2)*alpha[i]
        
        
        b,a = sig.butter(order, [W1, W2], btype='band')
        
        #Rescale decimated data
        if deci_rat > 1:
            decdata = sig.resample(data, int(len(data)/deci_rat))
        else:
            decdata = data
            
        placeholder = sig.lfilter(b,a,decdata)
    
        #Interpolate back up to original length of data
        #This ensures that the output filt_data is a nxm array, where
        #n is the number of center frequencies and m is the original length of 
        #the data
        if len(decdata) != len(data):
            dummy_time_act = np.arange(len(data))/fs
            dummy_time = np.arange(len(decdata))/fsdec
            placeholder = np.interp(dummy_time_act, dummy_time, placeholder)         
        
       
        #This initializes the filt_data array
        if i == 0:   
            filt_data = np.zeros((len(mid_bands),len(placeholder)))
        
        #Fill in filt_data with the filtered data held in placeholder
        for j in range(len(placeholder)):
            filt_data[i,j] = placeholder[j]
    
    return filt_data, mid_bands


def T60est(d,c = 1478,zi= 3.26e6,ai=0,alpha_p=0):
    """
    Parameters
    ----------
    d:          float;
                depth of water
    c:          float, Optional;
                speed of sound in water. Defaults to 1478m/s (rounded to nearest 
                whole value for any depth of water the tank can have), using 
                Garrett's eq. for the speed of sound in water relative to 
                temperature, depth, and salinity for a temparature of 19 degrees C 
                (rough avg. in tank).
    zi:         float, Optional;
                Acoustic impedance of side walls. Defaults to acoustic impdedance 
                of acrylic, accepted as 3.26E6 Ns/m**3 from the following source:
                https://www.ndt.net/links/proper.htm 
    ai:         float or ndarray of float, Optional;
                Absorption coefficient of tank walls. Defaults to 0 which ignores
                this input. If the absorption coefficient of the walls is known, 
                user can input this value and zi will be ignored, solving T60 
                using the known absorption. This may also be beneficial when 
                accounting for wall anechoic paneling (floor still assumed zi input). 
    alpha_p:    float or ndarray of float, Optional;
                Absorption coefficient due to thermoviscous molecular propagation
                losses. Defaults as 0 such that there is no propagation absorption.
                Can use alpha_prop(f,T,S,pH,depth) code to feed in an array of 
                frequency dependent absorption coefficients due to propagation 
                losses through the water.
                
    Returns
    -------
    T60:            float;
                    Estimate of the reverberation time (T60) in seconds. 
                    i.e. time it takes for the signal to drop by 60dB
    sigL:           float;
                    minimum excitation signal length (s) required by T60 based on 
                    Gemba recommendation for 5-10x length of T60. This gives 10x. 
    fschroeder:     float;
                    Schroeder Frequency (Hz). The lowest frequency of interest 
                    in which the tank is large.
    

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
    impedance is altered from the standard tank. May also input a wall absorption 
    coefficient if known. Such as from a measured wall absorption coefficient. 
    This currently still assumes the floor absorption is according to zi however.
    This improves all estimations given in this function.  
    
    Can also add in propagation absorption coefficients determined through 
    alpha_prop(f,T,S,pH,depth). Or leave that out by allowing the default to 
    remain 0. This further improves all estimations given in this function.
     
    last modified 9/1/2021      
    """
    import numpy as np
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    V = Lx*Ly*d #volume relative to current water depth
    A_floor = Lx*Ly #total surface area of tank floor
    A_acrylic = A_floor + 2*Ly*d + 2*Lx*d #total surface area of acrylic boundaries
    A_waterair = Lx*Ly #total surface area of water-air boundary
    S = A_acrylic +A_waterair #total surface area of (semi)absorptive boundaries
    
    #estimate absorption coefficients for boundaries
    zw = 1.5E6 #accepted acoustic impedance of water in Ns/m**3
    za = 415 #accepted acoustic impedance of air in Ns/m**3
    alpha_acrylic = 1-np.abs((zw-zi)/(zw+zi))
    alpha_air = 1-np.abs((zw-za)/(zw+za)) 
    Aw = alpha_air*A_waterair/S #water absorption coefficient spatially averaged
    
    if ai == 0:    
        #using zi (estimated acoustic impedance of walls)
        #Sum of alpha*A/S found in eq. 3.4.1 of Gemba
        #Absorption can be more thoroughly estimated using Physcs 661 notes. 
        Ai = alpha_acrylic*(A_acrylic)/S #acrylic absorp coeff spatially averaged
        Absorb = Aw + Ai #spatially averaged Absorption Coefficient
    else:
        #using ai (estimated acoustic absorption coefficient of walls) 
        Awall = ai*(A_acrylic-A_floor)/S
        Ai = alpha_acrylic*(A_floor)/S #acrylic absorp coeff spatially averaged
        Absorb = Aw + Ai + Awall #spatially averaged Absorption Coefficient
        
    #Eyring equation (661 notes eq. 4-24.124(reduce to 4-2.4.116 when alpha_p=0))
    T60 = (24*np.log(10)/c) * (V/(8*alpha_p*V - S*np.log(1-Absorb)))
    fschroeder = np.sqrt(c**3*T60/(V*4*np.log(10))) #Pierce eq6.6.4
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
        tstop = ((20+R)/60*texp*np.log(fstop/fstart) \
                 - l *np.log(fstop/f))/np.log(f/fstart)
    
    tls = l
    print(sigL)
    print(tstop)
    print(tls)
    return tstop, tls





def alpha_prop(f,T=16,S=5,pH=7.7,depth=0.6):
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
            depth of water in km. Defaults to 0.6m or 0.0006 km. Which will 
            make that term in the function negligible as basically zero. 
            Effective for 0<z<7000m or 0<z<7 km. 
    
    Returns
    -------
    a_p:    ndarray of float;
            absorption coefficient if sound (alpha) for propagation 
            losses through the water. (Np/m or Nepers/m)
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Primarily due to viscous effects above 100kHz (high), but also due to 
    Chemical relaxation of Boric Acid up to a few kHz (low), Chemical 
    relaxation of Magnesium Sulfate up to a few 100kHz (mid). This formulation 
    comes from Francois & Garrison 1982 Part1 and Part2 published in JASA
    Vol 72 No.3 and No.6 December 1982. See Specifically Fig. 7 Part 2
    
    Can apply this in propagation models similar to account for thermoviscous
    molecular losses.

    For reference: http://resource.npl.co.uk/acoustics/techguides/seaabsorption/
    
    last modified 9/1/2021      
    """
    import numpy as np
    f = f/1000 #convert freq input to kH
    c = 1412+3.21*T+1.19*S+0.0167*depth #m/s
    #boric acid term
    f1 = 2.8*(S/35)**0.5 *10**(4-1245/(T+273.1)) #kHz
    A1 = 8.86/c*10**(0.7*pH-5)
    P1 = 1 #pressure correction for boric acid not found important
    term1 = A1*P1*f1*f**2/(f1**2+f**2)
    #magnesium sulfate term
    f2 = 8.17*10**(8-1990/(T+273.1)) #kHz
    A2 = 21.44*S/c*(1+0.025*T) #dB/km/kHz
    P2 = 1 - 1.37e-4*depth + 6.2e-9*depth**2
    term2 = A2*P2*f2*f**2/(f2**2 + f**2)
    #pure water term
    if T <= 20:
        A3 = 4.937e-4 - 2.59e-5*T + 9.11e-7*T**2 - 1.5e-8*T**3 #dB/km/kHz**2
    else: 
        A3 = 3.964e-4 - 1.146e-5*T + 1.45e-7*T**2 - 6.5e-10*T**3 #dB/km/kHz**2
    P3 = 1 - 3.83e-5*depth + 4.9e-10*depth**2    
    term3 = A3*P3*f**2 #pure water term
    #print(f'Boron term={term1} in dB/km')
    #print(f'Magnesium term={term2} in dB/km')
    #print(f'Pure Water term={term3} in dB/km')
    a_p = (term1 + term2 + term3) 
    #Original function returns solution in dB/km
    #convert dB/km to Np/m (Nepers/meter)
    a_p = a_p/1000 *0.115129254650564 
    return a_p


#before the following function can be used, currently the variable below must
#be initialized to ensure it will function due to a conditional statement used.  
l1 = None #DO NOT ERASE
def T60meas_bounds(data,fs):
    """
    Parameters
    ----------
    data:       Ndarray;
                Impulse Response data.
    fs:         Float;
                Sampling rate. 
    
    Returns
    -------
    tbounds:    List (2 values);
                List of the initial and final time bounds to perform the 
                reverse Schroeder integration on T60meas(data,fs,t0,t1,d,c,rt,plot)

    Notes
    -----
    Author: Cameron Vongsawad
    
    Utilize a pop up graph to view the 10log10(h(t)**2) of an input impulse 
    response or h(t). This allows you to choose the appropriate time bounds 
    to perform the reverse Schroeder integration on the impulse response h(t)
    according to ISO354:2003, ISO3382-1:2009, ISO3382-2:2008 standards. 
    
    *****Because of the conditional statement in updatePlot(), the statement: 
    "l1=None" must remain before this function. This simply initializes l1 
    until another workaround is determined. 
    
    This is to be passed into the T60meas code in order to plot the decay curve 
    and determine the T60 of the impulse response. OctFilter() is recommended 
    prior to this function in order to pass through specific frequency bands.
    When doing this, it is also recommended that you write a loop to loop 
    through each octave band through this function.
    
    last modified 7/19/2021      
    """
    
    
    #Creatre fonts to be used in the plotting.
    LARGE_FONT= ("Verdana", 12)
    #Medium_FONT= ("Verdana", 10)
    
    #import necessary packages for use in GUI and plotting
    import tkinter
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
    # Implement the default Matplotlib key bindings.
    from matplotlib.backend_bases import key_press_handler
    from matplotlib.figure import Figure
    
    #Setup the main window for the GUI    
    root = tkinter.Tk()
    root.wm_title("Check h(t)**2")
    label = tkinter.Label(root, text="Impulse Response Squared in dB", 
                          font=LARGE_FONT)
    label.pack(pady=10,padx=10)
    #Create a figure to plot
    import numpy as np
    fig = Figure(figsize=(10, 8), dpi=100)
    ax = fig.add_subplot(111)#.plot(t,data)
    #Create time array for the sample data (IR)
    t = np.linspace(0,len(data)/fs,len(data))
    #Grab data to plot and adjust to properly view h(t)**2
    floor = np.max(data)*1e-5 
    ht = np.clip(data,floor,(np.max(data)+floor)) #data w/out zeros
    b = 10*np.log10((np.abs(ht))**2)
    #Plot initial data
    fig.suptitle("Choose a Time Interval (t1 to t2) for Reverse Schroeder" 
                 +"Integration. This should be where the plot is most linear.")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Level (dB)")
    ax.plot(t,b)
    ax.grid(True)
    fig.canvas.draw_idle()
    
    #Create inputs for choosing bounds on the graph
    in1_label = tkinter.Label(root,text="t1",font=LARGE_FONT)
    root.in1 = tkinter.Entry(root)
    in2_label = tkinter.Label(root,text="t2",font=LARGE_FONT)
    root.in2 = tkinter.Entry(root)
    in1_label.pack(side="left", fill="x", expand = False)
    in2_label.pack(side="right", fill="x", expand = False)
    root.in1.pack(side="left", fill="x", expand = False)
    root.in2.pack(side="right", fill="x", expand = False)
    
    #Replace vlines for checking bounds on the 10log10(h(t)**2) plot
    #when the Check button is pressed. 
    def updatePlot(t1,t2):
        tt1 = float(t1)
        tt2 = float(t2)
        global l1
        global l2 
        if l1 != None:
            ax.lines.remove(l1)
            ax.lines.remove(l2)
            fig.canvas.draw_idle()
        l1 = ax.axvline(tt1,color="orange",linestyle="--")
        l2 = ax.axvline(tt2,color="orange",linestyle="--")
        fig.canvas.draw_idle()
    
    #calculate the reverberation time based on the time bounds chosen when 
    #Run T60mean button is pressed.     
    def Reverb(t1,t2):
        tt1 = float(t1)
        tt2 = float(t2)
        global tbounds
        global t60
        tbounds = [tt1,tt2]
        root.destroy()
        
    #A tk.DrawingArea.
    canvas = FigureCanvasTkAgg(fig, master=root)  
    canvas.draw()
    
    #Create a toolbar for the GUI including ability to save plot. 
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas.mpl_connect("key_press_event", lambda event: print(
            f"you pressed {event.key}"))
    canvas.mpl_connect("key_press_event", key_press_handler)
    
    #Create buttons for control in GUI
    run_button = tkinter.Button(root,text="Run T60meas",command=lambda: Reverb(
            root.in1.get(),root.in2.get()))
    check_button = tkinter.Button(root,text="Check",command=lambda: updatePlot(
            root.in1.get(),root.in2.get()))
    
    # Packing order for Widgets are processed sequentially.
    # The canvas is rather flexible in its size, so we pack it last which makes
    # sure the UI controls are displayed as long as possible.
    run_button.pack(side="bottom")
    check_button.pack(side="bottom")
    toolbar.pack(side=tkinter.BOTTOM, fill="y")#tkinter.X)
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    #Loop GUI
    tkinter.mainloop()
    
    #Return the bounds in a list and the T60 time. 
    return tbounds

def T60meas(ht,fs,t0,t1,d=0.6,c=1478,rt='T60',plot=False):
    """
    Calculate the T20, T30, or T60 from Backward Schroeder Integration on the 
    measured impulse response hsys.
    
    Parameters
    ----------
    ht:         ndarray of float;
                Measured Impulse Response of the environment. 
    fs:         float;
                Sampling frequency of the impulse reponse. 
    t0:         int;
                start time in seconds
    t1:         int;
                finish time in seconds 
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
    

    Returns
    -------
    T60:        float;
                Calculated reverberation time (T60) in the tank in seconds.
                This is calculated using the Through The System (TTS) 
                response to evaluate reverberation only in the tank. 
                (i.e. time it takes for the signal in the tank to drop by 
                60dB)
   
    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the measured T60 in the tank. 
    
    Some guidance for this part found here: 
    https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/room.py
    the above link provides an alternate method to more generalize this solution 
    
    This also follows ISO3382-1:2009(E)
    
    last modified 5/18/2021      
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pylab
    from scipy import stats
    params = {'legend.fontsize': 24,
              'figure.figsize': (15, 10),
              'axes.labelsize': 28,
              'axes.titlesize':29,
              'axes.titleweight':'bold',
              'xtick.labelsize':24,
              'ytick.labelsize':24,
              'lines.linewidth':3}
    pylab.rcParams.update(params)
    
    ##avoid log(0) to prevent exploding by clipping all zeros to 0.00001% of 
    ##the max and go just beyond the max value so as to not clip. 
    floor = np.max(ht)*1e-5 
    ht1 = np.clip(ht,floor,(np.max(ht)+floor)) #data w/out zeros
    
    #Portion of array to actually look at for the T60meas. This is found by
    #eyeing it (ISO3382-1:2009(E)). Look just after the 10*np.log10(np.abs(ht)**2) 
    #is flat in the beg. and just before it is flat in the end (background level).
    t0 = int(t0*fs) #convert to samples
    t1 = int(t1*fs) #convert to samples
    ht1 = ht1[t0:t1] 

    #Backward Schroeder Integration
    T = 1/fs
    schroeder = np.cumsum(ht1[::-1]**2)[::-1]*T
    schroeder_dB = 10*np.log10(schroeder)  

    if rt == 'T10':
        #determine T10 between -5dB and -15dB of the max value of the decay curve
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
        #determine T60 between -5dB and -65dB of the max value of the decay curve
        init = -5.0
        end = -65.0
        factor = 1.0 #amount to mult. T60 by to extrapolate T60
    
    #Relative value to refine search for init & end bounds for rt measurement.
    maxval = np.max(schroeder_dB)
    schroeder_dB = schroeder_dB - maxval  
    #Linear regression
    #determine the value on the decay curve where it is nearest the init and end
    #values below the maximum of the the decay curve
    sch_init = schroeder_dB[np.abs(schroeder_dB - init).argmin()]
    sch_end = schroeder_dB[np.abs(schroeder_dB - end).argmin()]
    
    check_actual = (sch_init - sch_end +5)
    check_bounds = (init-end)
    if check_actual < check_bounds:
        raise ValueError(f"Decay not large enough for {rt} measurement."
                         +"Choose smaller rt value.")
        
    #indices of where the decay curve matches the init and end condition
    init_sample = np.where(schroeder_dB == sch_init)[0][0]
    end_sample = np.where(schroeder_dB == sch_end)[0][0]
       
    #Reverberation time (RT)
    #convert samples to time and determine the difference
    t_init = init_sample / fs
    t_end = end_sample / fs
    RT = t_end - t_init
    T60 = factor*RT
    print('T60 =',T60,'s') 
    print('')
    print('')
    
    
    if plot == True:
        t = np.linspace(0,len(ht1)/fs,len(ht1))
        Level = 10*np.log10((np.abs(ht1))**2)
        plt.figure()
        #plot the IR**2 in dB
        plt.plot(t,Level)
        plt.xlabel('Time (s)')
        plt.ylabel('Level (dB)')
        plt.grid()
        #plot Decay Curve
        plt.plot(t,(schroeder_dB + maxval))
        plt.legend([r'$10log[h^{2}(t)]$','Decay Curve'])
        est,_,_ = T60est(d,c)
        plt.title(f'T60meas={np.around(T60*1000,decimals=2)}ms,' 
                             +f'T60est={np.around(est*1000,decimals=2)}ms')    
    
    return T60



def alpha_wall(T60,d=0.6,c=1478,acc=False,alpha_p=0):
    """
    Calculate the spatially averaged absorption coefficient of the walls of the
    tank based on the measured T60 of the tank (either averaged or over freq)
    
    Parameters
    ----------
    T60:        float or array of float;
                Calculated reverberation time (T60) in the tank in seconds.
                This is calculated using the T60meas function.
    d:          float, Optional;
                depth of water. Defaults to a common 0.6m of water in the tank
    c:          float, Optional;
                speed of sound in water. Defaults to 1478 rounded to nearest whole
                value for any depth of water the tank can have, using Garrett's eq.
                for the speed of sound in water relative to temperature, depth, and 
                salinity for a temparature of 19 degrees C (rough avg. in tank).
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
    alpha_S:    float;
                Estimated spatially averaged absorption coefficient 
                (Nepers/m**2) for the room based on the measured T60 and the 
                Norris-Eyring Equation.  
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    Calculate the spatially averaged absorption coefficient of the tank boudaries
    from the measured T60 in the tank. 
    
    This assumes that the incident energy per unit time and area is the same for 
    all wall surfaces at any given time. If one of the consituent areas S1 has 
    an absorption coefficient alpha1 that is substantially different than the 
    coefficient alpha0 of the remaining surface area S-S1, this assumption 
    becomes questionable. See 661 notes section 4-2.4.2.4 (Rooms with 
    asymmetric or nonuniform absorption).
    
    
    last modified 10/21/2021      
    """
    import numpy as np
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    V = Lx*Ly*d #volume relative to current water depth
    S = 2*(Lx*Ly+Ly*d+d*Lx) #total enclosed surface area including air-water
    
    if acc == True:
        #account for perfectly reflective water-air boundary
        print('Calc. absorp. assuming pressure-release water-air (reflective)')
        Aw = 0
    else: 
        print('Calc. absorption w/ respect to small surface absorption')
        #account for slight impedance water-air boundary
        #estimate absorption coefficients for boundaries
        zw = 1.5E6 #accepted acoustic impedance of water in Ns/m**3
        za = 415 #accepted acoustic impedance of air in Ns/m**3
        alpha_air = 1-np.abs((zw-za)/(zw+za)) 
        Aw = alpha_air*Lx*Ly/(S) #water absorption coefficient spatially averaged    
        
    #solving Eyring equation w/propagation loss found in 661 notes eq 4-2.4.124
    #with default of alpha_p=0 this simplifies to eq. 4-2.4.116.
    alpha_S = 1-np.exp(V/S*(8*alpha_p - (24*np.log(10))/(c*T60)))
    
    #Account for the absorption of the water-air surface boundary while 
    #assuming that the incident energy per unit time and area is the same for 
    #all wall surfaces at any given time. (This does not apply asymmetric correction)
    A = alpha_S*S - Aw
    alpha_S = A/S
    return alpha_S 

def alpha_addition(ai,T60_1,T60_2,dS,d,c=1478):
    """
    Determine the absorption coefficient of any added material in the tank by
    comparison of a pre and post T60 measurement. This is performed by solving 
    for the change in absorptive area A between two measurements using the 
    Sabine equation. (Solving this for the Eyring Eq. may be a better idea)
    
    Parameters
    ----------
    ai:         float;
                Measured or estimated absorption of tank acrylic walls that are
                being covered by input material
    T60_1:      float;
                Measured reverbeation time in seconds of the tank prior to 
                application of new material into the tank. 
    T60_2:      float;
                Measured reverbeation time in seconds of the tank post  
                application of new material into the tank.
    dS:         float;
                Effective change in surface area (m**2) of the tank boundaries 
                due to application of new material. 
    d:          float;
                Depth of the water in m
    c:          float, Optional;
                Speed of sound in m/s. Defaults to 1478 m/s
    
    Returns:
    alpha_add:  float;
                Measured absorption (Nepers) of material input into the tank 
                based on reverberation time prior to placing new material.
        
    """
    import numpy as np
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    V = Lx*Ly*d #volume relative to current water depth
    #661 eq 4-2.4.90 & 4-2.4.91 generalized for any sound speed
    alpha_add = ai+24*np.log(10)*V/(c*dS)*(1/T60_2 - 1/T60_1) 
    return alpha_add

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
    #rigid wall solution for natural frequencies & 
    #Pressure release surface(nz component)
    fN = lambda nx,ny,nz: c/(2)*np.sqrt((nx/Lx)**2 + (ny/Ly)**2 + (
            (2*nz-1)/(2*Lz))**2)
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
    #order all the frequencies & associated modes in numerical order of freq.
    idxs = np.argsort(f) 
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
    #Eigen-Function for rectangular tank assuming rigid walls and pressure release 
    #water-air interface
    Psi = lambda nx,ny,nz :np.cos(nx*np.pi*x/Lx) * np.cos(ny*np.pi*y/Ly) * np.cos(
            (2*nz-1)*np.pi*z/(2*Lz))
    
    psi = []
    print('')
    print('calculating EigenFunctions')
    for i in range(len(mode)):
            psi.append(Psi(mode[i,0],mode[i,1],mode[i,2]))
    psi = np.array(psi)

    
    if len(x) > 1:
        ##########################################
        #The rest of this function is solely for #
        #contour plotting the Eigenfunctions psi #
        ##########################################
        if plot == True:
            print(f'plotting first {num} EigenFunctions')
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 24,
                      'figure.figsize': (15, 10),
                     'axes.labelsize': 28,
                     'axes.titlesize':29,
                     'axes.titleweight':'bold',
                     'xtick.labelsize':24,
                     'ytick.labelsize':24,
                     'lines.linewidth':3}
            pylab.rcParams.update(params)
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
                #check if mode actually present in this plane, if not, do not plot
                check =np.ones((len(x),len(y)))
                if np.any(psi1 != check) == True:
                    fig,ax=plt.subplots(1,1)
                    if pstyle == 'line':
                        cb = ax.contour(x,y,psi1,colors='black',linestyles='dashed')
                        ax.clabel(cb,inline=True,fontsize=15)
                    else: 
                        cb = ax.contourf(x,y,psi1)
                        fig.colorbar(cb)
                    ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)} Hz'
                                    + f'where Z={z}m')
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
                #check if mode actually present in this plane, if not, do not plot
                check =np.ones((len(x),len(z)))
                if np.any(psi2 != check) == True:
                    fig,ax=plt.subplots(1,1)
                    if pstyle == 'line':
                        cb = ax.contour(x,z,psi2,colors='black',linestyles='dashed')
                        ax.clabel(cb,inline=True,fontsize=15)
                    else:
                        cb = ax.contourf(x,z,psi2)
                        fig.colorbar(cb)
                    ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)} Hz' 
                                    + f'where Y={y}m')
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
                #check if mode actually present in this plane, if not, do not plot
                check =np.ones((len(y),len(z)))
                if np.any(psi3 != check) == True:
                    fig,ax=plt.subplots(1,1)
                    if pstyle == 'line':
                        cb = ax.contour(y,z,psi3,colors='black',linestyles='dashed')
                        ax.clabel(cb,inline=True,fontsize=15)
                    else:    
                        cb = ax.contourf(y,z,psi3)
                        fig.colorbar(cb)
                    ax.set_title(f'{modeint[i,:]} Mode f={np.round(fint[i],2)}' 
                                    + f'Hz where X={x}m')
                    ax.set_xlabel('Y (m)')
                    ax.set_ylabel('Z (m)')
                    plt.show()
                else: 
                    print('undesired mode not plotted in y-z')
    else:
        #not sure what this plot physically means, but I have it here for now. 
        #might need to change default to not plotting if singular value is input.
        if plot == True:
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 24,
                      'figure.figsize': (15, 10),
                     'axes.labelsize': 28,
                     'axes.titlesize':29,
                     'axes.titleweight':'bold',
                     'xtick.labelsize':24,
                     'ytick.labelsize':24,
                     'lines.linewidth':3}
            pylab.rcParams.update(params)
            plt.figure()
            plt.plot(f,psi)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel(r'$\Psi (r)$')
            plt.title(rf'Eigenfunction $\Psi$({x},{y},{z})')
    
    return psi 





def P_model_Pierce(Psi0,Psi,k,kn,mode,alpha,d=0.6,A=1,acc=False,anech=False,
                   alpha_p=0):
    """
    Parameters
    ----------
    Psi0:       Ndarray of float;
                Eigenfunctions at source position. Determined w/ TankFunc function.
    Psi:        Ndarray of float;
                Eigenfunctions at receiver position. Determined w/ TankFunc function.
    k:          Ndarray of float;
                Wavenumber of frequency band of interest.
    kn:         Ndarray of float;
                Eigenmodes of the tank environment. Determined w/ TankMode function
    mode:       Ndarray of float;
                Ordered Eigenmode array calculated and output by TankMode function
    alpha:      Ndarray of float;
                Spatially averaged absorption coefficient
    d:          float, Optional;
                depth of water. Defaults to a common 0.6m of water in the tank
    A:          float, Optional;
                Amplitude of the function. A=1 (default) ensures the solution is 
                simply the Green's function. 
    acc:        boolian, Optional;
                Account for the assumption that the water-air boundary is perfectly
                reflective. Defaults to False to not make this assumption and give
                the overall spatially averaged absorption coefficient. If True, 
                then the spatially averaged absorption coefficient that is returned
                only accounts for the walls and the floor of the enclosure. 
    anech:      boolian, Optional;
                Defaults to False to calculate with no anechoic panels in the tank. 
                If true, the calculation takes into account the thickness of the 
                panels on the inner dimensions of the tank environment for the 
                calculation of the spatially averaged absorption coefficient.
    alpha_p:    float or ndarray of float, Optional;
                Absorption coefficient due to thermoviscous molecular propagation
                losses. Defaults as 0 such that there is no propagation absorption.
                Can use alpha_prop(f,T,S,pH,depth) code to feed in an array of 
                frequency dependent absorption coefficients due to propagation 
                losses through the water. 
    
    Returns
    -------
    P:          Ndarray of float;
                Green's function or pressure in the tank as a function of source 
                and receiver positions. 
    
    Notes
    -----
    Author: Cameron Vongsawad
    
    This follows Pierce "Acoustics: An Introduction to its Physical Principles 
    and Applications" 3rd Ed. eq 6.5.20 (2019) or more precisely Leishman 661 
    notes 4-2C eq. 4-2.2.193 (2021)
    
    last modified 5/19/2021      
    """
    import numpy as np
    #dimensions of tank
    Lx = 1.22   #width of tank (m) 
    Ly = 3.66   #length of tank (m)
    #Alter the dimensions relative to the thickness of the anechoic panels
    if anech == True:
        Lx = Lx - 2*0.05
        Ly = Ly - 2*0.05
        print('Calculating w/ respect to anechoic panels')
    V = Lx*Ly*d #volume relative to current water depth
    if acc == True:
        S = 2*(Ly*d+d*Lx)+Lx*Ly #total enclosed surface area minus air-water surface
        print('Calculating w/ respect to walls only')
    else: 
        S = 2*(Lx*Ly+Ly*d+d*Lx) #total enclosed surface area including air-water
        print('Calculating w/ respect to water surface & walls')
    #Spatially averaged absorption area including propagation absorption 
    alpha_wall = alpha #wall absorption/impedance accounted for
    As = S*alpha_wall + 8*alpha_p*V 
    x,y,z = mode[0],mode[1],mode[2]
    if x == 0:
        Ex = 1
    else:
        Ex = 2
    if y == 0:
        Ey = 1
    else:
        Ey = 2
    Ez = 2
    lamb = 1/(Ex*Ey*Ez)
    P= -4*np.pi*A*np.sum((Psi*Psi0)/(V*lamb*(k**2-kn**2-1j*k*(As/(4*V)))))
    
    return P


