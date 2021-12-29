def freqVsRelTLwithAutospec(Gxx,f,start_freq = 0,end_freq = 0):
    '''
    This is a function that will calculate Transmission Loss relative to the first column
    in the Gxx array.

    Parameters :
    -----------
    x : ndarray of float
        Time series signal. Should be a single dimensional array. Input dtype
        should be float32 or float64.
    Gxx : ndarray of float
        A list of Gxx arrays returned from autospec on each position in the scan.
        It is and n x m array where n is the length of one Gxx list returned from
        autospec and m is the number of positions in the scan. If m = 1 this code
        will not work.
    f : ndarray of float
        The freq array returned by autospec. This is an n x 1 array.
    start_freq : int
        Most likely you only care about a certain band of frequencies, so here you
        can input the starting frequency. This will default to 0.
    end_freq : int
        Ending frequency for desired frequency band. Defaults to 0 but that will 
        just allow f to remain it's full length.

    Returns : 
    ---------
    rel_Gxx_level : ndarray of floats
        The relative Transmission loss array in dB. n x m-1 array
    freq_band : ndarray of floats 
        the f array that is cut off at the input frequencies
    '''
    if len(Gxx.shape) != 2:  # Check to see if Gxx is a 2-d array
        print('Gxx must be a 2 dimmensional array')
        return 
    if len(Gxx[:,0]) != len(f): # Check to see if they built the array correctly
        print('You built your Gxx array incorrectly or you did not use autospec')
        return 

    import numpy as np

    rel_Gxx_level = np.zeros((len(Gxx[:,0]),len(Gxx[0,:]) - 1)) # making empty array
    for i in range(len(Gxx[0]) - 1):
        rel_Gxx_level[:,i] = 10*np.log10(Gxx[:,i+1]/Gxx[:,0])  # filling it with the relative TL relative to the first position

    def findFreqIndex(f_array,freq): # function that finds the closest index to a desired frequency
        adjust_array = abs(f_array - freq)
        val = min(adjust_array)
        index = np.where(adjust_array == val)[0][0]
        return index

    if end_freq > f[-1]:  # checking to see if the ending frequency is in array
        print("You need a higher sampling frequency to include this ending frequency.")
        return

    if end_freq == 0: # if they want the freq range to be as high as possible then none of the frequency array needs to be chopped off at the end
        start_index = findFreqIndex(f,start_freq)
        end_index = len(f) - 1
    else:  # finding indices of the frequencies
        start_index, end_index = findFreqIndex(f,start_freq), findFreqIndex(f,end_freq)
    
    rel_Gxx_level = rel_Gxx_level[start_index:end_index,:] # reducing arrays to include what is desired. 
    freq_band = f[start_index:end_index]

    return rel_Gxx_level, freq_band

'''
# Using my function
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')
from autospec import autospec
from ESAUdata import ESAUdata
path = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-24B/2021-06-24_scan'
desire = [i for i in range(11)]
channels = [1]
fs = 1.5e6
sig_length = 0.65
N = int(fs*sig_length)
_,_,_,_,ch1,_,_ = ESAUdata(path,desire,channels,N,N)
ns = 2**15
unitflag = 0
pref = 1e-6
_,f,_ = autospec(ch1, fs, ns, N, unitflag, pref)
Gxx = np.zeros((ns//2,len(desire)))
for i in desire:
    Gxx[:,i],_,_ = autospec(ch1[:,i], fs, ns, N, unitflag, pref)

rel_Gxx_level, freq_band = freqVsRelTLwithAutospec(Gxx,f,10000,100000)

for i in desire[:-1]:
    plt.figure()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('TL relative to 0.1 meters')
    plt.plot(freq_band,rel_Gxx_level[:,i],label = str(0.1*i+0.2) + ' m')
    plt.legend()
    plt.show()
'''