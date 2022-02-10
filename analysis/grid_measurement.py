import numpy as np
import matplotlib.pyplot as plt
from ESAUdata import ESAUdata
from readLogFile import readLogFile
from ESAUpose import ESAUpose

path_sin = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-10-21/2021-10-21_scan'
# ALL HAVE 500 kHz SAMPLING
# Initial Scan:
# This is a 100 kHz scan, first half, 370 points, Use ESAUpose.

# Scan 1:
# This is an 80 kHz scan, first half, 314 points, Use ESAUpose.

# Scan 2:
# Second half of scan 1, 108 points, there is overlap, Use ESAUpose!

# Scan 3:
# Second half of Initial Scan, 27 points, there are some missing still,  Use ESAUpose.

# Scan4:
# 50 kHz, 405 points

# Scan 5:
# This one is 30 kHz, 405 points

# Scan 6:
# This one is 10 kHz, 405 points, weak signal (impedance matching transformer)

# Scan 7:
# This is 10 kHz, 405 points, strong signal (without Impedance Matching Transformer)

num_scan_sin = '' # the scan number in form of a string
num_points_sin = 0 # integer number of positions. NOTE a lot of them are different because we need to stitch multiple scans together

path_chirp = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-09-27/2021-09-27_scan'
# Initial Scan:
# BAD

# Scan 1:
# BAD

# Scan 2:
# BAD

# Scan 3:
# 5-10 kHz, 729 points

# Scan 4:
# 10k-50k, 729 points

# Scan 5:
# BAD, but if you really care: 50k-100k, 578 points, maybe some of it could be used if you use ESAUpose to ensure what positions there are.

# Scan 6:
# 50k-100k, 729 points

num_scan_chirp = '' # the scan number in form of a string. NOTE the bad ones
num_points_chirp = 729 # integer number of positions.

path = path_chirp + num_scan_chirp
num_points = num_points_chirp

desire = [i for i in range(len(num_points))]
channels=[0,1]
_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path)
N = int(fs*measurement_duration)
_,_,_,_,rec_sig,_,_ = ESAUdata(path, desire, channels, N, N)