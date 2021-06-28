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
    if len(Gxx.shape) != 2:
        return 'Gxx must be a 2 dimmensional array'
    if len(Gxx[:,0]) != len(f):
        return 'You built your Gxx array incorrectly or you did not use autospec'

    import numpy as np

    rel_Gxx_level = np.zeros((len(Gxx[:,0]),len(Gxx[0,:]) - 1))
    for i in range(len(Gxx[0,:]) - 1):
        rel_Gxx_level[:,i] = 10*np.log10(Gxx[:,i+1]/Gxx[:,0])

    def findFreqIndex(f_array,freq):
        adjust_array = abs(f_array - freq)
        val = min(adjust_array)
        index = np.where(adjust_array == val)[0][0]
        return index

    if end_freq == 0:
        start_index = findFreqIndex(f,start_freq)
        end_index = len(f) - 1
    else:
        start_index, end_index = findFreqIndex(f,start_freq), findFreqIndex(f,end_freq)
    
    rel_Gxx_level = rel_Gxx_level[start_index:end_index,:]
    freq_band = f[start_index:end_index]

    return rel_Gxx_level, freq_band