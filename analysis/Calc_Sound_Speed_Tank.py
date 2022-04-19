#%%
"""
Created on Wed May 26 16:47:56 2021

@author: scott
"""

def calc_sound_speed_tank(fs,rec_signal,signal_length_plus_zeros,rec_y_position,tolerance = 20.0):
    '''
    Code Description:
        This code takes a measurement and calculates the speed of sound 
        in the tank by comparing the delay between to received signals.
        This will not work if the measurement had only one position. Also,
        the receiver can only move in the y direction, and the source and
        receiver must be in the same point on the xz plane. 
        
    Parameters
    ----------
    fs : Float
        The sampling frequency of the measurement in Hz. This can be 
        obtained from the readLogFile function.
    rec_signal : ndarray
        The received signal from the measurement. This can be obtained from
        the ESAUdata funciton. This is array has dimensions (N,number of 
        scans). This function will return "Not able to calculate speed 
        of sound" if it is an N x 1 array.
    signal_length_plus_zeros : float
        This is the signal Duration plus the leading and trailing zeros. Signal
        Duration can be obtained from readLogFile function.
    rec_y_position : ndarray
        The y postions of the hydrophone at each scan. This can be obtained
        the readLogFile function, but you must construct the array. The array
        has dimensions (1,number of scans). 
    tolerance : float
        The tolerance level in error of the speed of sound. The function will 
        get rid of any values that have error greater than the tolerance and 
        will take the average of the remaining values. Defaults to 20 m/s
    

    Returns
    -------
    c : Float
        The speed of sound in the tank

    '''
    import numpy as np
    from scipy import signal
    import matplotlib.pyplot as plt
    if len(rec_signal[0,:])==1: 
        return print("Not able to calculate speed of sound")
        # checking to see if there are multiple scans
    if len(rec_y_position)!=len(rec_signal[0,:]):
        return print("Hydrophone positioning does not match number of scans")
        # checking to see if the positioning array was constructed correctly
    
    
    d = [] # this will be a list of distances between the positions
    for i in range(len(rec_y_position)-1):
        d.append(abs(rec_y_position[i+1]-rec_y_position[0]))
    
    delta_t = [] # this will be a list of all the times
    delta_t_error_min = []
    delta_t_error_max = []
    T = np.linspace(0,signal_length_plus_zeros,num = len(rec_signal[:,0]))
    Tcor = np.linspace(-signal_length_plus_zeros,signal_length_plus_zeros, num = 2*len(rec_signal[:,0])-1)
    for i in range(len(rec_y_position)-1):
        correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
        delta_t.append(Tcor[np.argmax(correlation)])
        delta_t_error_min.append(Tcor[np.argmax(correlation)-1])
        delta_t_error_max.append(Tcor[np.argmax(correlation)+1])
    
    c = [] # this will be a list of all the different calculated sound speeds
    c_error_min = []
    c_error_max = []
    for i in range(len(d)):
        c.append(d[i]/delta_t[i])
        c_error_min.append(d[i]/delta_t[i] - d[i]/delta_t_error_min[i])
        c_error_max.append(d[i]/delta_t_error_max[i] - d[i]/delta_t[i])
        
    c_minus_outlier = []
    c_error_min_minus_outlier = []
    c_error_max_minus_outlier = []
    d_minus_outlier = []
    for i in range(len(c)):
        if c[i] > 1000:
            c_minus_outlier.append(c[i])
            c_error_min_minus_outlier.append(c_error_min[i])
            c_error_max_minus_outlier.append(c_error_max[i])
            d_minus_outlier.append(d[i])
    
    c_error = [] # total error
    for i in range(len(c_error_min_minus_outlier)):
        c_error.append(c_error_min_minus_outlier[i]+c_error_max_minus_outlier[i])
    
    tolerated_error_index = []
    for i in range(len(c_error)):
        if c_error[i]<tolerance:
            tolerated_error_index.append(True)
        else:
            tolerated_error_index.append(False)
            
    d_tolerated = []
    c_tolerated = []
    c_error_min_tolerated = []
    c_error_max_tolerated = []
    d_excluded = []
    c_excluded = []
    c_error_min_excluded = []
    c_error_max_excluded = []
    
    for i in range(len(tolerated_error_index)):
        if tolerated_error_index[i]:
            d_tolerated.append(d_minus_outlier[i])
            c_tolerated.append(c_minus_outlier[i])
            c_error_min_tolerated.append(c_error_min_minus_outlier[i])
            c_error_max_tolerated.append(c_error_max_minus_outlier[i])
        else:
            d_excluded.append(d_minus_outlier[i])
            c_excluded.append(c_minus_outlier[i])
            c_error_min_excluded.append(c_error_min_minus_outlier[i])
            c_error_max_excluded.append(c_error_max_minus_outlier[i])
        
    plt.figure()
    plt.errorbar(d_tolerated,\
                 c_tolerated,\
                 yerr = [c_error_min_tolerated,c_error_max_tolerated],\
                 fmt = 'o',\
                 ecolor = 'green')
    plt.errorbar(d_excluded,\
                 c_excluded,\
                 yerr = [c_error_min_excluded,c_error_max_excluded],\
                 fmt = 'x',\
                 ecolor = 'red')
    plt.axhline(1478,linestyle='--',color = 'r')
    #plt.xlim(0.8,1.0)
    #plt.ylim(1430,1510)
    #plt.xlim(0,0.2)
    #plt.ylim(1430,1510)
    plt.show()
    
    return np.mean(c_minus_outlier)



import numpy as np
from ESAUdata import ESAUdata
from ESAUpose import ESAUpose
from readLogFile import readLogFile

path = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-06B/2021-08-06_scan1' # 16 points
num_points = 16
_,_,_,fs,_,_,_,measurementDuration,_,_,_,_,_,_,_ = readLogFile(f'/ID000_001log.txt',path)

'''
yRec = np.zeros(num_points)
for i in range(2):
    if i==0:
        for j in range(10):
            _,_,_,fs,_,_,_,measurementDuration,_,_,_,_,_,yRec[i*10+j],_ = readLogFile(f'/ID0{i}{j}_001log.txt',path)
    if i==1:
        for j in range(6):
            _,_,_,fs,_,_,_,measurementDuration,_,_,_,_,_,yRec[i*10+j],_ = readLogFile(f'/ID0{i}{j}_001log.txt',path)
'''

desire = [i for i in range(num_points)]
channels = [1]
_,_,dd = ESAUpose(path,desire)
_,_,_,_,rec_signal,_,_ = ESAUdata(path,desire,channels,fs*measurementDuration)

c = calc_sound_speed_tank(fs, rec_signal, measurementDuration, dd)

print(c)
# %%
