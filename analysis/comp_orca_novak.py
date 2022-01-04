import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import pickle

import sys


sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/uw-library')

import uwlib
from uwlib.orca import ORCA

#def forceAspect(ax,aspect=1):
#    im = ax.get_images()
#    extent =  im[0].get_extent()
#    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

SHOW_PLOTS = False
SAVE_PLOTS = True

logger = logging.getLogger(__name__)

SVP_FOLDER = "orcafiles"
SAVE_FOLDER = "defaultoutput/2021-10"

def main():

    max_modes = 1000

    OPT_FILE = 'orcafiles/orca_tank_opt.toml'
    test_env_files = ["svp_tank_air_0.5m.toml"] # hey make sure sound speed is right based on temperature!!!!!

    for testfile in test_env_files:
        full_file = os.path.join(SVP_FOLDER, testfile)

        logger.info(f"Telling ORCA to load {full_file}")
        orca = ORCA(base_svp=full_file, base_opt = OPT_FILE)

        # set transmission loss source depth (in m)
        src_depth = [0.25]   
        # set receiver depth(s) in m
        rec_depth = [0.25]
        # set the frequencies in Hz
        freqs = np.linspace(10000.0,50000.0,80001) # df of 0.5 Hz

        # set depths at which mode functions should be defined
        mode_depth = np.append(src_depth, rec_depth)

        # set ranges in m
        ranges = [0.3]
        #ranges_mirror = np.append(-1*np.flip(ranges),ranges)

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

    ORCA_dict = {'frequency':freqs,'transmission loss':tl[0,0,0,:]}
    filename = 'orca_tl_src&rec_depth_0.25m'
    outfile = open(filename,'wb')
    pickle.dump(ORCA_dict,outfile)
    outfile.close()
    
if __name__ == "__main__":

    logging.basicConfig(
        format='%(levelname)-8s %(asctime)s %(name)-8s %(module)s:%(funcName)s:%(lineno)d: %(message)s', level=logging.DEBUG)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    main()