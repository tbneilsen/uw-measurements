import logging
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

import pdb

import sys

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/uw-library')

def findFreqIndex(f_array,freq): # function that finds the closest index to a desired frequency
        adjust_array = abs(f_array - freq)
        val = min(adjust_array)
        index = np.where(adjust_array == val)[0][0]
        return index

def relOAPSLfft(Gxx,freq_index,plusorminus):
    listlen = len(Gxx[0,:])
    OASPL = np.zeros(listlen)
    for i in range(listlen):
        OASPL[i] = 10*np.log10(np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,i])/np.sum(Gxx[freq_index-plusorminus:freq_index+plusorminus,0]))
    return OASPL

import uwlib
from uwlib.orca import ORCA
from ESAUdata import ESAUdata
from readLogFile import readLogFile
from ESAUpose import ESAUpose
import TimeGate_UnderwaterTank as TG

import sys
#sys.path.append('/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20')
sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')

from autospec import autospec
from freq_vs_rel_TL_using_autospec import freqVsRelTLwithAutospec

path3 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-20_scan2' # 100 kHz, 100 points
path4 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-20_scan3' # 50 kHz, 100 points
path5 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20/2021-05-25_scan1' # 100 kHz, 100 points, Ran is on whole time

path6 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-07-09/2021-07-09_scan' # 100 kHz, 100 points, signal on longer
path7 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-07-09/2021-07-09_scan1' # 100 kHz, 100 points, signal on longer, longer settling time

path8 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-24/2021-06-24_scan7' # 1 - 10 kHz linear chirp, 11 points 
path9 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-06-24B/2021-06-24_scan' # 10 - 100 kHz linear chirp, 11 points

path10 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan' # 100 kHz
path11 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan1' # 80 kHz
path12 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan2' # 50 kHz
path13 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan3' # 30 kHz
path14 = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-03/2021-08-03_scan4' # 10 kHz

path_used = path14
freqs = [10000.0]
desire = [i for i in range(100)]
channel = [0,1]
c = 1478

_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path_used)
A, R, dd = ESAUpose(path_used,desire)
'''
acoustic_center_offset = 0.0092                 #This is not commented out for with_acoust_offset
dd = dd + 2*acoustic_center_offset              #This is not commented out for with_acoust_offset
for i in desire:                                #This is not commented out for with_acoust_offset
    A[i] = list(A[i])                           #This is not commented out for with_acoust_offset
    R[i] = list(R[i])                           #This is not commented out for with_acoust_offset
    A[i][1] = A[i][1] + acoustic_center_offset  #This is not commented out for with_acoust_offset
    R[i][1] = R[i][1] - acoustic_center_offset  #This is not commented out for with_acoust_offset
    A[i] = tuple(A[i])                          #This is not commented out for with_acoust_offset
    R[i] = tuple(R[i])                          #This is not commented out for with_acoust_offset
'''
tshort = np.zeros(len(desire))
tside = np.zeros(len(desire))
tdirect = np.zeros(len(desire))
for i in desire:
    tshort[i],tside[i],tdirect[i],_ = TG.gateValue(A[i], R[i], depth, c, 'tank',False)

N = int(fs*measurement_duration)

_, _, _, _, rec_sig, _, _ = ESAUdata(path_used,desire,channel,N)

gated_rec_sig = np.zeros((N,len(desire)))
for i in desire:
    gated_rec_sig[:,i] = TG.gatefunc(rec_sig[:,i],fs,tside[i],leading,0.1)


ns = 2**15
unitflag = 0
pref = 1e-6
Gxx = np.zeros((ns//2,len(desire)))
#OASPL = np.zeros(len(desire))
_, f, _ = autospec(rec_sig, fs, ns, N, unitflag, pref)
for i in desire:
    Gxx[:,i], _, _ = autospec(rec_sig[:,i], fs, ns, N, unitflag, pref)
'''
start_freq = 10000
end_freq = 100000
rel_Gxx_level, freq_band = freqVsRelTLwithAutospec(Gxx,f,start_freq,end_freq)


from Relative_TL_Tank import calcRelativeTransmissionLoss

rec_start = fs*(tdirect + leading)
rec_end = fs*(tside + leading)

rec_start = rec_start.astype('int')
rec_end = rec_end.astype('int')

rel_TL = calcRelativeTransmissionLoss(gated_rec_sig,rec_start,rec_end)
'''

plusorminus = 50
freq_index = findFreqIndex(f,freqs[0])
rel_TL_fft = relOAPSLfft(Gxx,freq_index,plusorminus)

SHOW_PLOTS = False
SAVE_PLOTS = True

logger = logging.getLogger(__name__)

SVP_FOLDER = "orcafiles"
SAVE_FOLDER = "/home/byu.local/sh747/underwater/scott-hollingsworth/orca_tank/2021-09"

# Plot TL for differenet frequencies vs which parameter?

plot_relative_tl = True # this is relative to the first point in the array
freq_vs_TL = False
range_vs_TL = True
depth_vs_TL = False


def main():

    if freq_vs_TL and range_vs_TL:
        print('You accidentally set freq_vs_TL and range_vs_tl both to TRUE. FIX IT!')
        return
    if freq_vs_TL and depth_vs_TL:
        print('You accidentally set freq_vs_TL and depth_vs_tl both to TRUE. FIX IT!')
        return
    if depth_vs_TL and range_vs_TL:
        print('You accidentally set depth_vs_TL and range_vs_tl both to TRUE. FIX IT!')
        return

    max_modes = 1000

    OPT_FILE = '/home/byu.local/sh747/underwater/scott-hollingsworth/underwater-measurements/analysis/orcafiles/orca_tank_opt.toml'
    test_env_files = ["svp_tank_air_0.519m.toml"]



    for testfile in test_env_files:
        full_file = os.path.join(SVP_FOLDER, testfile)

        logger.info(f"Telling ORCA to load {full_file}")
        orca = ORCA(base_svp=full_file , base_opt = OPT_FILE)
        # pdb.set_trace()

        # set transmission loss source depth (in m)
        src_depth = [depth/2]   # depth - A[0][2] used instead of depth/2 for chg_depth

        # set receiver depth(s) in m
        rec_depth = [depth/2]   # depth - R[18][2] used instead of depth/2 for chg_depth 

        # set depths at which mode functions should be defined
        # orca gets confused if the exact same number is in both src_depth and rec_depth
        mode_depth = np.append(src_depth, rec_depth)

        # set ranges in m
        ranges = dd #np.linspace(0.1, 1.1, 11) # [9.0]

        # set the source depth for the tl calculation
        logger.info("ORCA is set up and ready to go")

        # max modes for orca uses negative values to limit number of modes, but this does weird things
        # when trying to compare back to the rmin stuff
        # pdb.set_trace()
        # orca.max_modes(max_modes)
        #orca.set_num_modes = max_modes
        #orca.max_num_modes = max_modes

        # set depths at which mode functions should be defined
        orca.set_mode_grid(mode_depth, force_list=True)
        # run ORCA to get
        # kn = complex modal eignevalues, kn.shape = freq x mode number
        # phi = depth-dependent mode functions, phi,shape = freq x depth x mode number
        # nmode = the number of modes found at each frequency
        # mfz = depths at which mode functions are defined

        kn, phi, nmode, mfz = orca.run(freqs, parallel=False)

        zs_idx = []
        zr_idx =[] 

        # pdb.set_trace()
        # then calculate the transmission loss
        tl = uwlib.tl.calc_tl_from_orca(
            orca, freqs, ranges, src_depth, rec_depth,
            )

        if plot_relative_tl:
            if freq_vs_TL:
                for iir, r in list(enumerate(ranges))[1:]: # I am skipping the first range because that's what it's relative to
                    plt.figure()
                    plt.plot(freqs, tl[0,0,iir,:] - tl[0,0,0,:]*np.ones(len(tl[0,0,0,:])),\
                         label = "calc TL re " + str(ranges[0]) + " m @ " + str(r) + " m") # Calc Relative TL
                    # this is where you can plot data over the calculated transmission loss
                    # if you do plot data make sure you change save_name
                    plt.plot(freq_band,rel_Gxx_level[:,iir - 1], label = 'true TL re ' + str(ranges[0]) + ' m @ ' + str(r) + ' m')
                    plt.xlabel("Frequency, Hz")
                    plt.ylabel("TL, dB re 0.1 m") # change re depending on what you're doing
                    plt.title('Frequency vs Relative Transmission Loss\
                        \nReciever Depth: ' + str(rec_depth[0]) + ' m' +\
                        '\nRange: ' + str(r) + ' m' +\
                        '\n' + testfile[:-5]) 
                    plt.grid()
                    plt.legend()
                    plt.gcf().tight_layout()
                    # when creating the save name variable you must specify if you are plotting just ORCA, just data, or comparing the two
                    # orca_{etc} or data_{etc} or comp_{etc} 
                    # ALSO note the rtl instead of tl in the save_name!! When plot_rel_tl = False it is just 'tl'
                    save_name = 'comp_' + 'freq_vs_rtl_' + '@range' + str(r) + 'm' + '_' + testfile[:-5] + '.png'
                    plt.savefig(os.path.join( SAVE_FOLDER, save_name))
            if range_vs_TL:  #if plotting range vs. TL
                for iif, f in enumerate(freqs):
                    # Do you want all the freqs on the same graph? 
                    plt.figure()
                    # You also need to comment/uncomment out plt.savefig() in or out of loop!!
                    plt.plot(ranges, tl[0,0,:,iif] - tl[0,0,0,iif]*np.ones(len(tl[0,0,:,iif])),\
                         label = "calc TL re " + str(ranges[0]) + ' m @ ' + str(f) + " Hz") # Calc Relative TL
                    # this is where you can plot actual data over the calculated tl
                    plt.plot(ranges,rel_TL_fft, label = "true TL re " + str(ranges[0]) + ' m @ ' + str(f) + " Hz") # true Relative TL
                    plt.xlabel('Range, m')
                    plt.ylabel("TL, dB re "+str(round(dd[0],3))+" m") # change re depending on what it's to
                    plt.title('Range vs Relative Transmission Loss\
                        \nReciever Depth: ' + str(rec_depth[0]) + ' m' +\
                        '\nFrequency: ' + str(f) + ' Hz' +\
                        '\n' + testfile[:-5]) 
                    plt.grid()
                    plt.legend()
                    plt.gcf().tight_layout()
                    # If you want all freqs to plot on the same graph change @freq to a frequency band or just get rid of it
                    # when creating the save name variable you must specify if you are plotting just ORCA, just data, or comparing the two
                    # orca_{etc} or data_{etc} or comp_{etc} 
                    # ALSO note the rtl instead of tl in the save_name!! When plot_rel_tl = False it is just 'tl'
                    save_name = 'comp_' + 'range_vs_rtl_' + '@freq' + str(f) + 'Hz' + '_' + testfile[:-5] + '_sum+-' + str(plusorminus) + 'bins.png'
                    #plt.savefig(os.path.join( SAVE_FOLDER, save_name))
                    plt.savefig(os.path.join( SAVE_FOLDER, save_name))
            if depth_vs_TL: #if plotting rec_depth vs. TL
                print('unsure of which depth to calc TL relative to')
                # Maybe one day I will come in and add stuff here?
            
        else:
            if freq_vs_TL:
                for iir, r in enumerate(ranges):
                    plt.figure()
                    plt.plot(freqs, tl[0,0,iir,:],\
                         label = "calc TL " + str(r) + " m") # Calc Absolute TL
                    # this is where you can plot data over the calculated transmission loss
                    # if you do plot data make sure you change save_name
                    #plt.plot(freq_band,Gxx_level[:,iir - 1], label = 'true TL ' + str(r) + ' m')
                    plt.xlabel("Frequency, Hz")
                    plt.ylabel("TL, dB re 1 \mu Pa")
                    plt.title('Frequency vs Transmission Loss\
                        \nReciever Depth: ' + str(rec_depth[0]) + ' m' +\
                        '\nRange: ' + str(r) + ' m' +\
                            '\n' + testfile[:-5]) 
                    plt.grid()
                    plt.legend()
                    plt.gcf().tight_layout()
                    # when creating the save name variable you must specify if you are plotting just ORCA, just data, or comparing the two
                    # orca_{etc} or data_{etc} or comp_{etc} 
                    save_name = 'orca_' + 'freq_vs_tl_' + '@range' + str(r) + 'm' + '_' + testfile[:-5] + '.png'
                    plt.savefig(os.path.join( SAVE_FOLDER, save_name))
            if range_vs_TL:  #if plotting range vs. TL
                for iif, f in enumerate(freqs):
                    # Do you want all the freqs on the same graph? 
                    plt.figure()
                    # You also need to comment/uncomment out plt.savefig() in or out of loop!!
                    plt.plot(ranges, tl[0,0,:,iif],\
                         label = "calc TL " + str(f) + " Hz") # Calc TL
                    # this is where you can plot actual data over the calculated tl
                    #plt.plot(Range,TL, 'o', markersize = 1, label = "true TL " + str(f) + " Hz") # true TL
                    plt.xlabel('Range, m')
                    plt.ylabel("TL, dB re 1 \mu Pa")
                    plt.title('Range vs Transmission Loss\
                        \nReciever Depth: ' + str(rec_depth[0]) + ' m' +\
                        '\nFrequency: ' + str(f) + ' Hz' +\
                            '\n' + testfile[:-5])
                    plt.grid()
                    plt.legend()
                    plt.gcf().tight_layout()
                    # If you want all freqs to plot on the same graph change @freq to a frequency band or just get rid of it
                    # when creating the save name variable you must specify if you are plotting just ORCA, just data, or comparing the two
                    # orca_{etc} or data_{etc} or comp_{etc} 
                    save_name = 'orca_' + 'range_vs_tl_' + '@freq' + str(f) + 'Hz' + '_' + testfile[:-5] + '.png'
                    plt.savefig(os.path.join( SAVE_FOLDER, save_name))
                #plt.savefig(os.path.join( SAVE_FOLDER, save_name))
            if depth_vs_TL: #if plotting rec_depth vs. TL
                for iif,f in enumerate(freqs):
                    # Do you want all the freqs on the same graph? 
                    # plt.figure()
                    # You also need to comment/uncomment out plt.savefig() in or out of loop!!
                    plt.figure()
                    plt.plot(rec_depth, tl[0,:,0,iif], label=str(f) + " Hz")
                    plt.xlabel("Receiver Depth, m")
                    plt.ylabel('TL, dB re 1 \mu Pa')
                    plt.title('Depth vs Transmission Loss\
                            \nRange: ' + str(rec_depth[0]) + ' m' +\
                            '\nFrequency: ' + str(f) + ' Hz' +\
                                '\n' + testfile[:-5])
                    plt.grid()
                    plt.legend()
                    # If you want all freqs to plot on the same graph change @freq to a frequency band or just get rid of it
                    # when creating the save name variable you must specify if you are plotting just ORCA, just data, or comparing the two
                    # orca_{etc} or data_{etc} or comp_{etc} 
                    save_name = 'orca_' + 'range_vs_tl_' + '@freq' + str(f) + 'Hz' + '_' + testfile[:-5] + '.png'
                    save_name =testfile[:-5]+"_rel_tl_zr_"+str(ranges)+"m.png"
                    #plt.savefig()
                #plt.savefig()
                
                    
    #ORCA_dict = {'frequency':freqs,'transmission loss':tl,'ranges':ranges}
    #filename = 'ORCA_TL'
    #outfile = open(filename,'wb')
    #pickle.dump(ORCA_dict,outfile)
    #outfile.close()



if __name__ == "__main__":

    logging.basicConfig(
        format='%(levelname)-8s %(asctime)s %(name)-8s %(module)s:%(funcName)s:%(lineno)d: %(message)s', level=logging.DEBUG)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    main()
# %%
