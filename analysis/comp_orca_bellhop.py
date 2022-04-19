#import pickle
import matplotlib.pyplot as plt
import os
import sys
import logging
import numpy as np
# filename = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/bellhop/freq_tl_bellhop'
# infile = open(filename,'rb')
# bellhop_dict = pickle.load(infile)
# infile.close()
# freqs = bellhop_dict['frequency'] # these freqs will also be used in ORCA
# tl_freqs_bellhop = bellhop_dict['transmission loss']

sys.path.append('/home/byu.local/sh747/virtualenvs/bellhop/lib/python3.6/site-packages')
import arlpy.uwapm as pm

hw = 0.5
max_range = 1.5
c = 1484
zs = hw/2
zr = hw/2
r = 0.5
freqs = np.linspace(50000.0,100000.0,200)

# bathy = [r,hw]
# r : horizontal distance from transmitter
# hw : water depth
# gives the water depth at multiple points in the tank
bathy = [
    [0.0, hw],  # water depth at the transmitter
    [max_range, hw], # water depth 1 m away from transmitter
]

# ssp = sound speed profile [z,c]
# z = depth
# c = sound speed at given depth (z)
ssp = [
    [0.0, c],  # c m/s at the surface
    [hw/3, c],  # c m/s at hw/3 depth
    [hw/2, c],  # c m/s at hw/2 depth
    [hw, c],  # c m/s at hw depth (bottom)
]

# define environment

env = pm.create_env2d(
    depth=bathy,
    soundspeed=ssp,
    bottom_soundspeed=2750,
    bottom_density=1190, #
    bottom_absorption=0.128, # dB/mkHz
    tx_depth= zs,
    rx_depth= zr,
    rx_range= r, # m
    frequency = freqs
)

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/bellhop')
from plot_bellhop import calc_tl_freq

tl_freq_bellhop = calc_tl_freq(env)

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/uw-library')

import uwlib
from uwlib.orca import ORCA

SHOW_PLOTS = False
SAVE_PLOTS = True

logger = logging.getLogger(__name__)

SVP_FOLDER = "orcafiles"
SAVE_FOLDER = "defaultoutput/2022-02"

def main():

    max_modes = 1000
    hw = 0.5

    OPT_FILE = 'orcafiles/orca_tank_opt.toml'
    test_env_files = ["svp_tank_air_0.5m.toml"] # hey make sure sound speed is right based on temperature!!!!!

    for testfile in test_env_files:
        full_file = os.path.join(SVP_FOLDER, testfile)

        logger.info(f"Telling ORCA to load {full_file}")
        orca = ORCA(base_svp=full_file, base_opt = OPT_FILE)

        # set transmission loss source depth (in m)
        src_depth = zs
        # set receiver depth(s) in m
        rec_depth = zr

        # set depths at which mode functions should be defined
        mode_depth = np.append(src_depth, rec_depth)

        # set ranges in m
        ranges = r

        # set the source depth for the tl calculation
        logger.info("ORCA is set up and ready to go")

        # max modes for orca uses negative values to limit number of modes, but this does weird things
        # when trying to compare back to the rmin stuff

        # set depths at which mode functions should be defined
        orca.set_mode_grid(mode_depth, force_list=True)
        # run ORCA to get
        # kn = complex modal eignevalues, kn.shape = freq x mode number
        # phi = depth-dependent mode functions, phi,shape = freq x depth x mode number
        # nmode = the number of modes found at each frequency
        # mfz = depths at which mode functions are defined

        kn, phi, nmode, mfz = orca.run(freqs, parallel=False)

        zs_idx = [] # ASK DR. NEILSEN WHAT THESE ARE
        zr_idx = []

        # then calculate the transmission loss
        tl = uwlib.tl.calc_tl_from_orca(
            orca, freqs, ranges, src_depth, rec_depth,
            )
        
        # tl[src_depth, rec_depth, range, frequency]
        tl_freq_orca = tl[0,0,0,:]
        plt.figure()
        plt.plot(freqs/1000, tl_freq_bellhop)
        plt.plot(freqs/1000, tl_freq_orca)
        plt.legend(['bellhop','orca'])
        plt.xlabel('Frequency, kHz')
        plt.ylabel("Transmission Loss, dB")
        plt.grid()
        plt.title(
            'Comparing Orca To Bellhop: ' + str(round(freqs[0]/1000,2)) + '-' + str(round(freqs[-1]/1000,2)) + ' kHz' '\n'
            'Source Depth: ' + str(src_depth) + ' m' + '\n'
            'Receiver Depth: ' + str(rec_depth) + ' m' + '\n'
            'Source Receiver range: ' + str(ranges) + 'm'
            )
        plt.gcf().tight_layout()
        save_name =\
          'comp_orca_bellhop_'\
        + str(round(freqs[0]/1000,2)) + 'to' + str(round(freqs[-1]/1000,2)) + 'kHz_at_'\
        + 'range_' + str(ranges) + 'm_'\
        + 'src_depth_' + str(src_depth) + 'm_'\
        + 'rec_depth_' + str(rec_depth) + 'm_'\
        + testfile[:-5] + '.png'
        plt.savefig(os.path.join( SAVE_FOLDER, save_name))
        

if __name__ == "__main__":

    logging.basicConfig(
        format='%(levelname)-8s %(asctime)s %(name)-8s %(module)s:%(funcName)s:%(lineno)d: %(message)s', level=logging.DEBUG)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    main()