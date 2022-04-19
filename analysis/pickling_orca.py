import logging
import os
import numpy as np
import pickle
import sys

sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/uw-library')

import uwlib
from uwlib.orca import ORCA

logger = logging.getLogger(__name__)

SVP_FOLDER = "orcafiles"

max_modes = 1000

OPT_FILE = 'orcafiles/orca_tank_opt.toml'
test_env_files = ['svp_tank_air_0.47m.toml']


for testfile in test_env_files:
    water_level = testfile.lstrip('svp_tank_air_')
    water_level = water_level.rstrip('m.toml')
    water_level = float(water_level)
    
    full_file = os.path.join(SVP_FOLDER, testfile)

    logger.info(f"Telling ORCA to load {full_file}")
    orca = ORCA(base_svp=full_file, base_opt = OPT_FILE)
    # water depth
    hw = 0.47 # m
    # set source depth(s) in m
    src_depth = np.array([hw/2])
    # set receiver depth(s) in m
    rec_depth = np.array([hw/2])
    # set ranges in m
    ranges = np.linspace(0.1, 1.6, 150) # same as several measurements with sine waves including 2022-02-24
    # set the frequencies in Hz
    freqs = np.array([71000.0,82000.0,100000.0,115000.0,133000.0,200000.0])  # same as several measurements with sine waves including 2022-02-24
    # set depths at which mode functions should be defined
    mode_depth = np.append(src_depth, rec_depth)
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

    # then calculate the transmission loss
    tl = uwlib.tl.calc_tl_from_orca(
        orca, freqs, ranges, src_depth, rec_depth,
        )
    # tl[src_depth, rec_depth, ranges, freqs]
    orca_dict = {
        'src_depth':src_depth, 
        'rec_depth': rec_depth, 
        'ranges': ranges,
        'frequency':freqs,
        'transmission loss':tl
        }
    # tl[src_depth, rec_depth, ranges, freqs]
save_path = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2022-02/'
filename = 'modeling_paper_orca_pickle'
outfile = open(save_path+filename,'wb')
pickle.dump(orca_dict,outfile)
outfile.close()