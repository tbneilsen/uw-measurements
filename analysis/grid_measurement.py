#%%
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

num_scan_chirp = '6' # the scan number in form of a string. NOTE the bad ones
num_points_chirp = 729 # integer number of positions.

path = path_chirp + num_scan_chirp
num_points = num_points_chirp

desire = [i for i in range(num_points)]
channels=[0,1]
_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path)
N = int(fs*measurement_duration)
_,_,_,_,rec_sig,_,_ = ESAUdata(path, desire, channels, N, N)
plot = False
Acal=(0.6, 2.14, depth/2)
Rcal=(0.6, 2.06, depth/2)
A,R,dd = ESAUpose(path, desire, plot, Acal, Rcal)
#%%
chirp_grid = np.zeros((27,27))
for i in range(27):
    for j in range(27):
        chirp_grid[i,j] = 27*j + i
chirp_grid = chirp_grid.astype('int')
comp_array = [3,12,21]
for i in chirp_grid[3,comp_array]:
    print('Aegir: ' + str(A[i]))
for i in chirp_grid[3,comp_array]:
    print('Ran: ' + str(R[i]))
for i in chirp_grid[3,comp_array]:
    print('Range: ' + str(dd[i]))
## %%
import sys
sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')
from autospec import autospec
ns = 2**15
unitflag = 0
pref = 1e-6
Gxx = np.zeros((ns//2,3))
#OASPL = np.zeros(len(desire))
_, f, _ = autospec(rec_sig, fs, ns, N, unitflag, pref)
for i,ii in enumerate(chirp_grid[3,comp_array]):
    Gxx[:,i], _, _ = autospec(rec_sig[:,ii], fs, ns, N, unitflag, pref)

from freq_vs_rel_TL_using_autospec import freqVsRelTLwithAutospec
start_freq = 50000
end_freq = 100000
rel_Gxx_level, freq_band = freqVsRelTLwithAutospec(Gxx,f,start_freq,end_freq)
save_path = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2022-02/'
re = str(round(dd[chirp_grid[3,comp_array[0]]],2)) + 'm'
for i,ii in enumerate(chirp_grid[3,comp_array]):
    filename = 'frequency_at_' + str(round(dd[ii],2)) + 'm_re_' + re + '.jpg'
    plt.figure()
    plt.plot(freq_band,rel_Gxx_level[:,i])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('TL vs. Frequency at ' + str(round(dd[ii],2)) + '\n' + 're ' + re)
    plt.savefig(save_path + filename)
# %%
