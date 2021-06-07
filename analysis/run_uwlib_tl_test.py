import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
import matplotlib.gridspec as gridspec

import pdb

import sys

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/uw-library')

import uwlib
from uwlib.orca import ORCA
from ESAUdata import ESAUdata
from readLogFile import readLogFile

SHOW_PLOTS = False
SAVE_PLOTS = True

logger = logging.getLogger(__name__)

SVP_FOLDER = "orcafiles"
SAVE_FOLDER = "defaultoutput/2021-06"
# OPT_FILE  = "Hadassah.toml"

# Plot TL for differenet frequencies vs which parameter?
# plot_xaxis = "rec_depth"
plot_xaxis = "range"
plot_relative_tl = True # this is relative to the first point in the array
scale_factor = 1 # this is how much we are scaling up the tank


def main():

    max_modes = 1000

    # ["svp.gsl.no_shear.toml", "svp.gsl.toml"]
    test_env_files = ["svp_tank.toml"]

    freqs = [100000.0/scale_factor] #np.array([300.0, 250.0, 200.0, 100.0])

    for testfile in test_env_files:
        full_file = os.path.join(SVP_FOLDER, testfile)

        logger.info(f"Telling ORCA to load {full_file}")
        orca = ORCA(base_svp=full_file) # , base_opt = OPT_FILE)
        # pdb.set_trace()

        # set transmission loss source depth (in m)
        src_depth = [0.121*scale_factor]#np.linspace(0.05, 5.0, 512)  # 1.0

        # set receiver depth(s) in m
        rec_depth = [0.121*scale_factor] #[2.0,3.0, 3.5]#np.array([1.0, 2.0])

        # set depths at which mode functions should be defined
        # orca gets confused if the exact same number is in both src_depth and rec_depth
        mode_depth = np.append(src_depth, rec_depth)

        # set ranges in m
        ranges = np.linspace(0.1*scale_factor, 1.1*scale_factor, 1000) # [9.0]

        # set the source depth for the tl calculation
        logger.info("ORCA is set up and ready to go")

        # max modes for orca uses negative values to limit number of modes, but this does weird things
        # when trying to compare back to the rmin stuff
        # pdb.set_trace()
        orca.max_modes(max_modes)
        orca.set_num_modes = max_modes
        orca.max_num_modes = max_modes

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

        plt.figure()
        for iif, f in enumerate(freqs):
            # plot TL in order: src_depth, rec_depth, range, freq
            # select source depth, receiver depth, range, freq
            # pdb.set_trace()
            if plot_relative_tl:
                if plot_xaxis == "range": #if plotting TL vs. range
                    plt.plot(ranges, tl[0,0,:,iif] - tl[0,0,0,iif]*np.ones(len(tl[0,0,:,iif])), label=str(f) + " Hz") # Relative TL
                    plt.xlabel("Range, m")
                    plt.title(testfile[:-5] + " " + str(rec_depth[0]) + ' m' + f' scaled by {scale_factor}') 
                    save_name =testfile[:-5]+"_rel_tl_zr_"+str(rec_depth[0])+"m_" + str(f) + f"Hz_scaled_by_{scale_factor}.png" 
                elif plot_xaxis == "rec_depth": #if plotting TL vs. rec_depth
                    plt.plot(rec_depth, tl[0,:,0,iif], label=str(f) + " Hz")
                    plt.plot(rec_depth, tl[0,:,0,iif] - tl[0,0,0,iif]*np.ones(len(tl[0,0,:,iif])), label=str(f) + " Hz") # Relative TL
                    plt.xlabel("Receiver Depth, m")
                    plt.title(testfile[:-5] + " " + str(ranges[0]) + ' m') 
                    save_name =testfile[:-5]+"_rel_tl_zr_"+str(ranges)+"m.png" 
            else:
                if plot_xaxis == "range": #if plotting TL vs. range
                    plt.plot(ranges, tl[0,0,:,iif], label=str(f) + " Hz") # Absolute TL
                    plt.xlabel("Range, m")
                    plt.title(testfile[:-5] + " " + str(rec_depth[0]) + ' m') 
                    save_name =testfile[:-5]+"_tl_zr_"+str(rec_depth[0])+"m_" + str(f) + f"Hz_scaled_by_{scale_factor}.png" 
                elif plot_xaxis == "rec_depth": #if plotting TL vs. rec_depth
                    plt.plot(rec_depth, tl[0,:,0,iif], label=str(f) + " Hz") # Absolute TL
                    plt.xlabel("Receiver Depth, m")
                    plt.title(testfile[:-5] + " " + str(ranges[0]) + ' m') 
                    save_name =testfile[:-5]+"_tl_zr_"+str(ranges)+"m.png" 

        plt.ylabel("TL, dB re 1 \mu Pa")
        
        plt.grid()
        plt.legend()
        plt.gcf().tight_layout()
        plt.savefig(os.path.join( SAVE_FOLDER, save_name))
        print(kn[0,0].real*ranges[0])


if __name__ == "__main__":

    logging.basicConfig(
        format='%(levelname)-8s %(asctime)s %(name)-8s %(module)s:%(funcName)s:%(lineno)d: %(message)s', level=logging.DEBUG)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    main()
