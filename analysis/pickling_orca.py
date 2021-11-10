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
test_env_files = ['svp_tank_air_0.191m.toml','svp_tank_air_0.209m.toml','svp_tank_air_0.242m.toml','svp_tank_air_0.4m.toml','svp_tank_air_0.45m.toml','svp_tank_air_0.47m.toml','svp_tank_air_0.5m.toml','svp_tank_air_0.55m.toml','svp_tank_air_0.6m.toml','svp_tank_air_0.509m.toml','svp_tank_air_0.519m.toml']
dict_to_be_pickled = {}

for testfile in test_env_files:
    water_level = testfile.lstrip('svp_tank_air_')
    water_level = water_level.rstrip('m.toml')
    water_level = float(water_level)
    
    full_file = os.path.join(SVP_FOLDER, testfile)

    logger.info(f"Telling ORCA to load {full_file}")
    orca = ORCA(base_svp=full_file, base_opt = OPT_FILE)

    # set transmission loss source depth (in m)
    off_top = 0.01 # take 1 cm off the top and bottom for the source and receiver depths
    src_depth = np.linspace(off_top ,water_level - off_top, 481)   # this has 1 mm seperation if the depth is 0.5 meters. Also halfway is always included!
    # set receiver depth(s) in m
    rec_depth = np.linspace(off_top ,water_level - off_top, 481) # same as src_depth array
    # set ranges in m
    ranges = np.linspace(0.01, 3, 2991) # 10 cm to 3 meters with 1 mm seperation
    # set the frequencies in Hz
    freqs = np.linspace(5000.0,10000.0,5001)  # 1 kHz to 5 kHz with df of 1 Hz 
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
    dict_to_be_pickled.update({str(water_level)+' m':{'src_depth':src_depth,'rec_depth':rec_depth,'ranges':ranges,'freqs':freqs,'tl':tl}})
    # tl[src_depth, rec_depth, ranges, freqs]
  
filename = 'pickled_orca_1k_5k'
outfile = open(filename,'wb')
pickle.dump(dict_to_be_pickled,outfile)
outfile.close()