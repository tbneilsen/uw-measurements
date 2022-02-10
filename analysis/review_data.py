#%%
from ESAUdata import ESAUdata
from ESAUpose import ESAUpose
from readLogFile import readLogFile
from TimeGate_UnderwaterTank import uwsoundspeed
from pathlib import Path
import numpy as np

from run_uwlib_tl_test import SAVE_FOLDER
#import antigravity
'''
D = 1
T = 1
S = 0.9
model = 'Garrett'
c = uwsoundspeed(D,T,S,model)
'''
path = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-10-21/2021-10-26_scan7'

desire = [i for i in range(150)]
channels=[1]
_,_,_,fs,leading,signal_duration,trailing,measurement_duration,depth,_,_,_,_,_,_ = readLogFile('/ID000_001log.txt',path)

N=int(fs*measurement_duration)
_,_,cal,tegam,rec_sig,_,_ = ESAUdata(path,desire,channels,N,N)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(rec_sig[100000:100100,0])
#plt.show()
plt.figure()
plt.plot(rec_sig[100000:100100,50])
#plt.show()
plt.figure()
plt.plot(rec_sig[100000:100100,100])
#plt.show()
plt.figure()
plt.plot(rec_sig[100000:100100,149])
#plt.show()

'''
plot = False
Acal = (0.6,2.14,0.25)
Rcal = (0.6,2.06,0.25)
A, R, dd = ESAUpose(path,desire,plot,Acal,Rcal)

index_array = np.zeros((9,45)) # Aegir's Grid is 3x3 (9 points) and Ran's Grid is 3x15 (45 points)
for i in range(len(index_array[:,0])):
  for j in range(len(index_array[0,:])):
    index_array[i,j] = 9*j+i  # make a list of indices that jives with the way ESAU does scans. 
index_array=index_array.astype(int)

source_pos = []
for i in index_array[:,0]:
    source_pos.append(A[i][1:])

sr_dist = np.zeros((9,45))  # Aegir's Grid is 3x3 (9 points) and Ran's Grid is 3x15 (45 points)
for i in range(len(index_array[:,0])):
    for j,jj in enumerate(index_array[i,:]):
        sr_dist[i,j] = dd[jj]
# Do we want to make TL relative to calibration position?
# Or maybe to the first point in every grid?
# First point in every row?
# Try to go for Absolute TL?? Why is SL sometimes smaller than RL
import matplotlib.pyplot as plt
for i in range(4):
    plt.figure()
    plt.plot(cal[:,i])
'''
'''
import matplotlib.pyplot as plt
fig = plt.figure()
AEgir = fig.add_subplot(111)
Ran = fig.add_subplot(111)
lineA, = AEgir.plot(Acal[1],Acal[2],'.')
lineR, = Ran.plot(Rcal[1],Rcal[2],'.')
plt.ylim((0,0.5))
plt.xlim((0,3.66))
fig.show()
for i in range(50):
    lineA.set_xdata(A[i][1])
    lineA.set_ydata(A[i][2])
    lineR.set_xdata(R[i][1])
    lineR.set_ydata(R[i][2])
    fig.show()
'''
'''
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 24,
          'figure.figsize': (15, 10),
         'axes.labelsize': 28,
         'axes.titlesize':29,
         'axes.titleweight':'bold',
         'xtick.labelsize':24,
         'ytick.labelsize':24,
         'lines.linewidth':3}
pylab.rcParams.update(params)
import TimeGate_UnderwaterTank as TG
from ESAUpose import ESAUpose
import os
SAVE_FOLDER = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2021-11'
save_name = 'time_gate3'
c=1486.54900429149
A, R, dd=ESAUpose(path,desire,plot = False,Acal = (0.6,2.14,0.3), Rcal = (0.6,2.06,0.3))
tshort,tside,tdirect,_ = TG.gateValue(A[0], R[0], depth, c, 'tank', False)
gated_rec_sig = TG.gatefunc(rec_sig[:,0],fs,tside,leading,0.1)
plt.figure()
chop = rec_sig[int(0.05*fs):,0]
t = np.linspace(0,measurement_duration,N)
#plt.plot(t,rec_sig[:,0])
plt.plot(t,gated_rec_sig)
plt.axvline(leading+tdirect,linestyle='--',color = 'r')
plt.axvline(leading+tside,linestyle='--',color = 'r')
plt.title('Time Gating Example')
plt.xlabel('Time')
plt.ylim(-0.5,0.5)
plt.xlim(0.0475,0.0575)
plt.ylabel('Amplitude')
plt.savefig(os.path.join(SAVE_FOLDER, save_name))
'''

# %%
