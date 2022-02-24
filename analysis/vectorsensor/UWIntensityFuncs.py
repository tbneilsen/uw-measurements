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
from numpy.fft import fft



def intensityTimeAveraged(p, v, channel, fs, sig_cond = 0.01, rho = 1023, c = 1478):

    N = len(v)
    df = fs/N
    
    #Replace zeros with very low non-zero values, to avoid having NaN's in 
    #calculations
    v = np.where(v!=0,v,1e-16)
    p = np.where(p!=0,p,1e-16)/sig_cond
    
    v = fft(v)[0:N//2-1]*np.sqrt(2/N/df)
    
    #Take FFT of data to get time-averaged intensity as a function of frequency
    p = fft(p)[0:N//2 - 1]*np.sqrt(2/N/df)
    
    fss = np.linspace(0,(fs//2)-df,len(v))
    v = PMD_sensitivities(fss,v,channel)
    
    v = v/rho/c
    
    #Calculate intensity in dB re I_ref
    I = 1/2*np.conj(p)*v
    
    return I, fss


def PMD_sensitivities(fss,data,channel):
    """
    Applies the known sensitivities or frequency responses of the 
    particle motion detector. Outputs the adjusted pressure values for each
    component of the particle velocity in addition to the omnidirectional
    sensor channel. 
    Inputs should already be in frequency domain.
    Only gives accurate values between 10 Hz and 5 kHz.

    Parameters
    ----------
    fss : array of float64
        The single-sided frequency array over which we have data.
    data : array of float64
        Frequency response for particle velocity in x-direction,
        in Volts.
    channel : string
        Should be "x", "y", "z", or "omni". No other value will be accepted. 
        This indicates which sensitivity graph you want to use. 

    Returns
    -------
    data: array of float64
        Sensitivity-adjusted data in Pa. 
    

    Author: Corey Dobbs
    Sensitivity curves generated from data from GeoSpectrum, INC.
    """
    
    #Raw sensitivity data
    x_freq = [10.22029846850511, 11.232411940394085, 12.078660286617875, 13.083352347643194, 14.274928933739758, 15.575029292857643, 17.24221127581675, 19.227002329372812, 21.131047318035975, 23.392956474498433, 25.52348127773076, 27.84806008855136, 31.508129118298015, 34.37774718744459, 38.05760933395839, 42.43850928451851, 47.66869505038102, 53.1559494377222, 60.142266323764765, 67.06538251985371, 75.87978211352387, 84.61448036513036, 95.04249629341362, 109.9029432716683, 124.34755006595034, 140.6905418389522, 166.27241724092744, 199.38136798786553, 242.58167190591448, 278.4806008855136, 324.3702103077194, 369.6778234396754, 418.26455976002546, 473.2373104604571, 551.220032675185, 628.2133886344558, 705.635743329074, 781.1679579532339, 909.8937391768545, 1000, 1115.112326476353, 1234.4761584615364, 1386.6148584567668, 1626.8841634450591, 1854.124220183775, 2097.8128286144606, 2533.877552279372, 2339.2956474498424, 2951.425870949254, 3220.2284068670087, 3513.512263500637, 3833.507088952023, 4243.8556321224305, 4596.85316584366, 5052.07972687205, 4836.62655933848]
    x_sens = [-208.4070869653136, -206.78467022727608, -205.01475185167556, -203.39233511363804, -201.76991837560053, -200, -197.93511210765865, -196.01770334729017, -194.24780185088224, -192.4778778488842, -190.85546111084668, -189.2330387464116, -187.16814522767274, -185.5457228632377, -183.77581574043222, -181.85840698006376, -180.08849985725828, -178.3185871080553, -176.40118397408435, -174.9262576208148, -173.1563448716118, -171.3864321224088, -169.76401538437128, -167.9941026351683, -166.81416267783214, -165.48672670933053, -163.5693235753596, -161.65192044138863, -159.5870212962522, -158.11209494298265, -156.34218219377965, -155.01475185167556, -153.68731588317394, -152.50737592583778, -150.44247678070138, -148.52507364673045, -146.60766769956072, -144.9852509615232, -143.95280701535253, -144.8377577635565, -146.1651965452569, -147.49262407416222, -149.41003002133192, -151.62241673803751, -153.68731588317394, -155.60471901714487, -158.25958532775059, -156.9321549856465, -160.47197485765494, -161.7994108261566, -163.1268467946582, -164.30678675199437, -166.37168027073324, -169.3215329772724, -173.7463120370811, -171.53392250717675]

    y_freq = [9.909758078848219, 11.353196822953851, 12.430489852679903, 13.610006096712409, 15.037143481096349, 16.765233849708785, 18.691905945505994, 21.029783733210216, 23.660053676799187, 26.379101408912007, 29.145198758141483, 32.79051523862945, 37.56672670191928, 42.26533099460339, 47.98466215690239, 54.97402758126427, 65.30715445041264, 71.50408236224351, 81.91924709403268, 91.33346277987756, 100.91063697452304, 114.56583231014056, 128.89507810072928, 146.33711911875895, 167.65233849708804, 192.0722969866029, 220.0492256799066, 245.33743260978574, 278.5364160457071, 319.10756332960057, 365.5882358931654, 422.6533099460337, 493.0753414379461, 559.7981399890175, 635.5498464418795, 689.5769670386085, 761.8856098178795, 872.8605901390437, 1000, 1094.8889591650495, 1243.0489852679889, 1361.000609671241, 1531.2258027706503, 1676.522225442225, 1920.721641442376, 2141.4537621097566, 2431.2346048487034, 2810.7287163773094, 3162.2765666100563, 3655.8823589316503, 4076.0216693792845, 4382.602323984733, 4669.728150666063]
    y_sens = [-207.9044171528419, -205.51470572766962, -203.67647048511307, -201.83823524255652, -200, -197.79413033109188, -195.77206034458217, -193.5661766510521, -191.54411367685336, -189.52207173958746, -187.86764319174029, -185.6617665105212, -183.45588281699116, -181.43381984279242, -179.4117779055265, -177.3897079190168, -174.63235505518196, -173.5294097022615, -171.50735374037367, -170.03676694879064, -168.56617314489665, -166.91176563398233, -165.99264801270408, -164.52206122112102, -162.86764669789576, -161.76470835728622, -160.11029383406094, -158.8235277618092, -157.5367616895574, -156.06617489797438, -154.59558810639135, -153.30882203413958, -151.47058679158306, -149.81617928066873, -147.79411630646996, -146.32352951488693, -145.2205911742774, -144.3014735529991, -145.2205911742774, -146.5073502342182, -147.97794403811218, -149.81617928066873, -151.28676607225177, -152.94117358316606, -154.7794088257226, -156.61764406827916, -158.27205859150442, -160.47794228503443, -162.49999824692227, -164.3382334894788, -166.91176563398233, -169.48529076617493, -172.97794053195673]

    z_freq = [9.566862388340175, 10.569119902030973, 11.676377371085806, 13.043226578054817, 14.896258189732244, 16.825257726100954, 19.215613052315987, 21.94554634072577, 25.063317151616765, 28.624049154731658, 33.054518844441, 37.750556578990135, 43.59363963788101, 49.78696460168567, 60.09609595054467, 68.63391601717404, 76.66827553621107, 87.56042230152339, 101.11323027004049, 119.37774172866301, 137.8551927523018, 162.75656066360915, 187.94840367110714, 226.86591542855678, 270.82717785452724, 316.2277660168373, 373.34932805697815, 426.39086652181646, 503.4116564485246, 594.3450381934508, 663.9197032031976, 749.894367669458, 895.2062053979245, 1022.3868065816192, 1167.6367509419492, 1378.5519275230165, 1574.400729562301, 1778.2786590867531, 2076.3835157654607, 2293.9126173145755, 2619.8053378674313, 2959.0610385986183, 3455.105498823412, 3945.970480875676, 4311.3685045425345, 4816.061906793983]
    z_sens = [-207.13237153411274, -205.1470572766961, -202.94117498562827, -200.9558859725311, -198.52941040349256, -196.1029600787735, -193.89707778770565, -191.470593803894, -189.48531320556995, -186.83823594378765, -184.63235365271981, -182.20588649845442, -180.2205890705841, -177.79412191631872, -175.1470614840827, -173.16176405621238, -171.6176531842834, -169.63235575641306, -167.86765160651345, -166.10294745661378, -164.7794130331092, -163.01470888320955, -161.4705895965074, -159.48530058341026, -157.72058801873746, -155.95588386883782, -154.63235786010637, -153.08823857340423, -151.3235344235046, -149.11765213243675, -147.13235470456644, -145.3676505546668, -144.04412454593535, -145.58823541786433, -147.13235470456644, -149.77941513680247, -151.54411928670214, -153.30882343660176, -155.95588386883782, -157.2794182923424, -159.04412244224204, -161.02941145533924, -163.45588702437774, -165.66176931544558, -168.3088297476816, -172.94117919301485]

    omni_freq = [9.701568819841693, 10.732518029004055, 12.11526958519891, 13.95521945327884, 16.074603096500134, 18.515843945286544, 20.901396413198633, 24.075698173692196, 28.29790445436345, 32.59551625751754, 37.926900369292674, 45.487758256657514, 52.92782178361213, 65.43185547589711, 80.88999066075529, 92.23848046958804, 110.6266236760992, 135.3876262377098, 170.78746260371736, 211.13563346950744, 240.75698173692174, 294.6444099474015, 379.2690036929264, 459.4949765443655, 534.650841480977, 591.4661743862587, 654.3190589097071, 769.0683808259566, 876.9657766362001, 1051.791917044903, 1236.2478501095281, 1528.3075693434448, 1796.3301110585737, 2243.24786677543, 2583.9306380103762, 2916.840233091969, 3292.638844431824, 3716.85443613376, 3949.040752729703, 4153.569143802017, 4503.072819873462, 4784.369332225705]
    omni_sens = [-187.77777349641227, -186.36363075893968, -184.74747093316728, -182.72726922433728, -181.11112481148075, -179.29292477803475, -178.0808222482358, -176.86868889260512, -175.45454615513253, -174.24242050595979, -173.03030256324496, -171.6161598257724, -170.40404188305757, -169.19191623388483, -167.97979829117003, -167.17171452505485, -166.56565555369744, -166.16161367063987, -165.75757949404021, -165.55555469928245, -165.15151281622485, -164.74747863962523, -164.34343675656763, -163.93939487351003, -163.73737778521024, -163.93939487351003, -163.53535299045245, -163.53535299045245, -163.33333590215264, -163.93939487351003, -164.14141196180987, -164.14141196180987, -164.34343675656763, -164.74747863962523, -163.93939487351003, -162.72727693079526, -161.31313419332267, -159.4949495727925, -157.8787897470201, -156.06060512648995, -157.07070598090493, -158.8888906014351]
    
    #Interpolate up to frequecy length and apply sensitivities, based on
    #what input is given
    
    if channel == "x":
        x_sens = np.interp(fss,x_freq,x_sens)
        x_sens = 1e6*10**(x_sens/20)
        data = data/x_sens
    elif channel == "y":
        y_sens = np.interp(fss,y_freq,y_sens)
        y_sens = 1e6*10**(y_sens/20)
        data = data/y_sens
    elif channel == "z":
        z_sens = np.interp(fss,z_freq,z_sens)
        z_sens = 1e6*10**(z_sens/20)
        data = data/z_sens
    elif channel == "omni":
        omni_sens = np.interp(fss,omni_freq,omni_sens)
        omni_sens = 1e6*10**(omni_sens/20)
        data = data/omni_sens
    
    return data


def getIntensityBlocks(vx,vy,vz,p,ns):
    
    
    
    #Replace zeros with 1e-16
    vx = np.where(vx!=0, vx, 1e-16)
    vy = np.where(vy!=0, vy, 1e-16)
    vz = np.where(vz!=0, vz, 1e-16)
    p = np.where(p!=0, p, 1e-16)
    
    #Find the block FFT of each channel
    vxf,num_Blocks1 = cbFFT.computeBlockFFT(vx,ns)
    vyf,num_Blocks2 = cbFFT.computeBlockFFT(vy,ns)
    vzf,num_Blocks3 = cbFFT.computeBlockFFT(vz,ns)
    pf,num_Blocks4 = cbFFT.computeBlockFFT(p,ns)
    
    #Calculate the frequency dependent instantaneous intensity
    Ix = np.conj(pf)*vxf
    Iy = np.conj(pf)*vyf
    Iz = np.conj(pf)*vzf
    
        
    Ix_active = np.real(Ix)
    Ix_reactive = np.imag(Ix)
    Iy_active = np.real(Iy)
    Iy_reactive = np.imag(Iy)
    Iz_active = np.real(Iz)
    Iz_reactive = np.imag(Iz)
    
    #Put in a decibel scale
    I_ref = 6.61e-19 #Based on standard reference pressure in water, 1 microPa
    
    Ix_active = 10*np.log10(np.abs(Ix_active)/I_ref)
    Ix_reactive = 10*np.log10(np.abs(Ix_reactive)/I_ref)
    Iy_active = 10*np.log10(np.abs(Iy_reactive)/I_ref)
    Iy_reactive = 10*np.log10(np.abs(Iy_reactive)/I_ref)
    Iz_active = 10*np.log10(np.abs(Iz_active)/I_ref)
    Iz_reactive = 10*np.log10(np.abs(Iz_reactive)/I_ref)
    
    
    return Ix_active, Ix_reactive, Iy_active, Iy_reactive, Iz_active, Iz_reactive


def PMD_processing(fs,vx,vy,vz,omni,rho0=1023,c=1478):
    """
    NOTE: Obsolete.
    
    Use intensityTimeAveraged instead. 
    
    
    This code inputs the data from the 4 channels of the particle motion
    detector, and outputs the processed data. Code includes application of
    frequency-dependent sensitivities, and conversion from decibels to 
    the physical parameters of Pascals and m/s

    Parameters
    ----------
    fs : float
        The sampling frequency of the data.
    vx : array of float64
        The unprocessed measurements of particle velocity in the x-direction.
    vy : array of float64
        The unprocessed measurements of particle velocity in the y-direction.
    vz : array of float64
        The unprocessed measurements of particle velocity in the z-direction.
    omni : array of float64
        The unprocessed measurements of pressure from the omnidirectional
        sensor. 
    rho0 : float, optional
        The density of the medium. The default is 1023kg/m^3, for water.
    c : float, optional
        The speed of sound in the medium. The default is 1478m/s, for water.

    Returns
    -------
    vx : array of float64
        Time-domain x-direction velocity with adjustments for frequency response.
        In m/s
    vy : array of float64
        Time-domain y-direction velocity with adjustments for frequency response.
        In m/s
    vz : array of float64
        Time-domain z-direction velocity with adjustments for frequency response.
        In m/s
    p : array of float64
        Time-domain omni-directional pressure with adjustments for frequency response.
        In m/s
    
    Author: Corey Dobbs
    Last modified: 9/27/21
    References autospec.py in byuarglib.
    ***NEEDS TO BE UPDATED BASED ON PMD_sensitivities***

    """

    
    pref = 1e-6
    
    #Find frequency spectrum of data
    spec_x,fss_x,OASPLx = byu.autospec(vx,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
    spec_y,fss_y,OASPLy = byu.autospec(vy,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
    spec_z,fss_z,OASPLz = byu.autospec(vz,fs,ns=2**15,N=0,unitflag=0,pref=1e-6)
    spec_omni,fss_omni,OASPLp = byu.autospec(omni,fs,ns=2**15,N=0,unitflag=0,pref=1e-6) 
    
    #Convert to dB re 1 microPa or dB re 1 V (equivalent)
    vx_dB = 20*np.log10(spec_x/pref)
    vy_dB = 20*np.log10(spec_y/pref)
    vz_dB = 20*np.log10(spec_z/pref)
    omni_dB = 20*np.log10(spec_omni/pref)
    
    x_data, y_data, z_data, omni_data = PMD_sensitivities(fss_x,vx_dB,vy_dB,vz_dB,omni_dB)
    
    vx = pref*10**(x_data/20)/rho0/c #m/s
    vy = pref*10**(y_data/20)/rho0/c #m/s
    vz = pref*10**(z_data/20)/rho0/c #m/s
    p = pref*10**(omni_data/20) #Pa
    
    return vx, vy, vz, p



    
    
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
    

    
    
    
    