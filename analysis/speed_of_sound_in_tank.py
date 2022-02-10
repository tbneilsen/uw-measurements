#%%
import os
import numpy as np
from ESAUdata import ESAUdata
from ESAUpose import ESAUpose
from readLogFile import readLogFile

path = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-08-06B/2021-08-06_scan1' # 16 points
num_points = 3
_,_,_,fs,leading_zeros,_,trailing_zeros,measurementDuration,_,_,_,_,_,_,_ = readLogFile(f'/ID000_001log.txt',path)
tolerance = 20 # m/s
N = int(fs*measurementDuration)
desire = [i for i in range(num_points)]
channels = [1]
_,_,dd = ESAUpose(path,desire)
_,_,_,_,rec_signal,_,_ = ESAUdata(path,desire,channels,N)
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

d = [] # this will be a list of distances between the positions
for i in range(len(dd)-1):
    d.append(abs(dd[i+1]-dd[0]))

delta_t = [] # this will be a list of all the times
delta_t_error_min = []
delta_t_error_max = []
T = np.linspace(0,measurementDuration,num = len(rec_signal[:,0]))
Tcor = np.linspace(-measurementDuration,measurementDuration, num = 2*len(rec_signal[:,0])-1)
for i in range(len(dd)-1):
    correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
    delta_t.append(Tcor[np.argmax(correlation)])
    delta_t_error_min.append(Tcor[np.argmax(correlation)-1])
    delta_t_error_max.append(Tcor[np.argmax(correlation)+1])
    


# We're gonna hardcode i = 1 because correlate is not very good on this point
i = 1
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)]) #-390
plt.figure()
plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
plt.show()
delta_t[i] = shift/fs
'''
# We're gonna hardcode i = 2 because correlate is not very good on this point
i = 2
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-390
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 4 because correlate is not very good on this point
i = 4
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])+390
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 6 because correlate is not very good on this point
i = 6
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-1580
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 8 because correlate is not very good on this point
i = 8
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-1520
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 9 because correlate is not very good on this point
i = 9
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-3600
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 10 because correlate is not very good on this point
i = 10
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-2335
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 11 because correlate is not very good on this point
i = 11
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])+810
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 12 because correlate is not very good on this point
i = 12
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-2000
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 13 because correlate is not very good on this point
i = 13
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-1170
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

# We're gonna hardcode i = 14 because correlate is not very good on this point
i = 14
correlation = signal.correlate(rec_signal[:,i+1], rec_signal[:,0], mode="full")
shift = int(fs*Tcor[np.argmax(correlation)])-1525
#plt.figure()
#plt.plot(rec_signal[int(fs*leading_zeros)+5000+shift:N-int(fs*trailing_zeros)+shift,i+1])
#plt.plot(rec_signal[int(fs*leading_zeros)+5000:N-int(fs*trailing_zeros),0])
#plt.ylim(-0.1,0.1)
#plt.show()
delta_t[i] = shift/fs

c = [] # this will be a list of all the different calculated sound speeds
c_error_min = []
c_error_max = []
for i in range(len(d)):
    c.append(d[i]/delta_t[i])
    c_error_min.append(d[i]/delta_t[i] - d[i]/delta_t_error_min[i])
    c_error_max.append(d[i]/delta_t_error_max[i] - d[i]/delta_t[i])
    
c_minus_outlier = []
c_error_min_minus_outlier = []
c_error_max_minus_outlier = []
d_minus_outlier = []
for i in range(len(c)):
    if c[i] > 1000:
        c_minus_outlier.append(c[i])
        c_error_min_minus_outlier.append(c_error_min[i])
        c_error_max_minus_outlier.append(c_error_max[i])
        d_minus_outlier.append(d[i])

c_error = [] # total error
for i in range(len(c_error_min_minus_outlier)):
    c_error.append(c_error_min_minus_outlier[i]+c_error_max_minus_outlier[i])

tolerated_error_index = []
for i in range(len(c_error)):
    if c_error[i]<tolerance:
        tolerated_error_index.append(True)
    else:
        tolerated_error_index.append(False)
        
d_tolerated = []
c_tolerated = []
c_error_min_tolerated = []
c_error_max_tolerated = []
d_excluded = []
c_excluded = []
c_error_min_excluded = []
c_error_max_excluded = []

for i in range(len(tolerated_error_index)):
    if tolerated_error_index[i]:
        d_tolerated.append(d_minus_outlier[i])
        c_tolerated.append(c_minus_outlier[i])
        c_error_min_tolerated.append(c_error_min_minus_outlier[i])
        c_error_max_tolerated.append(c_error_max_minus_outlier[i])
    else:
        d_excluded.append(d_minus_outlier[i])
        c_excluded.append(c_minus_outlier[i])
        c_error_min_excluded.append(c_error_min_minus_outlier[i])
        c_error_max_excluded.append(c_error_max_minus_outlier[i])
    
plt.figure()
plt.errorbar(d_tolerated,\
                c_tolerated,\
                yerr = [c_error_min_tolerated,c_error_max_tolerated],\
                fmt = 'o',\
                ecolor = 'green',\
                label = 'measured sound speed')
plt.errorbar(d_excluded,\
                c_excluded,\
                yerr = [c_error_min_excluded,c_error_max_excluded],\
                fmt = 'x',\
                ecolor = 'red')
plt.axhline(1483.65,linestyle='--',color = 'r',label = 'predicted sound speed') # Garret model prediction
plt.title('Calculating the Speed of Sound in the Tank')
plt.xlabel('Range (m)')
plt.ylabel('Sound Speed (m/s)')
plt.legend()
#plt.xlim(0.8,1.0)
#plt.ylim(1430,1510)
#plt.xlim(0,0.2)
plt.ylim(1470,1515)
#plt.show()
SAVE_FOLDER = '/home/byu.local/sh747/underwater/scott-hollingsworth/codes/underwater-measurements/analysis/defaultoutput/2021-11'
save_name = 'speed_of_sound_in_tank_2.png'
plt.savefig(os.path.join( SAVE_FOLDER, save_name))
#print(np.mean(c_minus_outlier))
'''
# %%
