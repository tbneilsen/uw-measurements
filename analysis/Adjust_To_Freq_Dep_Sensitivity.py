def adjust(rec_sig,transducer,folder):
    '''
    Code Description:
        
    Parameters
    ----------
    rec_sig : Array of Floats
        This is the received signal from ESAU. It can be a 1 dimmensional array
        or a two dimmensional array that has multiple scans.   
    transducer : String
        Input which transducer was selected. For now it can only be a BK8103 or
        a TC4038
    folder : String
        Folder that contains all the Transducer Files

    Returns
    -------
    adjusted_rec_sig : Array of Floats
        Now rec_sig has been adjusted with the frequency dependant hydrophone
        sensitivities. 
    '''

    import numpy as np
    from scipy.interpolate import interp1d

    n = np.shape(rec_sig)[1]
    BK = False
    Teledyne = False
    filename = folder + transducer + '.txt'

    if transducer == 'BK8103-3249167':
        BK = True
    elif transducer == 'BK8103-3249178':
        BK = True
    elif transducer == 'BK8103-3249187':
        BK = True
    elif transducer == 'BK8103-3249189':
        BK = True
    elif transducer == 'BK8103-3249190':
        BK = True
    elif transducer == 'TC4038-0420006':
        Teledyne = True
    elif transducer == 'TC4038-0420008':
        Teledyne = True
    elif transducer == 'TC4038-4318006':
        Teledyne = True
    else:
        return print('Do not have hydrophone sensitivity data for the specified hydrophone. Check your spelling :)')
    
    f = open(filename,'r')
    lines = f.readlines()
    sens_string = lines[8].split('\t')
    freq_string = lines[11].split('\t')
    sensitivity_dB = np.zeros(len(sens_string))
    frequency = np.zeros(len(freq_string))
    for i in range(len(sens_string)):
        sensitivity_dB[i] = float(sens_string[i])
        frequency[i] = float(freq_string[i])
    sensitivity = (10**(sensitivity_dB/20))*(10**6)
    interpolation = interp1d(frequency,sensitivity)
    freq_min = frequency[0]
    freq_max = frequency[-1]

    for i in range(n):
        np.fft.fft(rec_sig[:,i])
