# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:21:13 2020

TEST CODE

@author: cvong
"""

import numpy as np
import byuarglib as byu
import SystemFResponse as sys
import MeasuredGreens as meas
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 15,
          'figure.figsize': (15, 10),
         'axes.labelsize': 24,
         'axes.titlesize':28,
         'axes.titleweight':'bold',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large',
         'lines.linewidth':2}
pylab.rcParams.update(params)
import scipy.signal as sci 
import ESAUpose as pose
import ESAUdata as data
import TimeGate_UnderwaterTank as tg
#import pdb
#pdb.set_trace()

###############################################################################
###############################################################################
#when switching file name via copy/paste, must change forward and back slashes
###############################################################################
###############################################################################
date =  '2021-02-11' #scan 14,15,16 2 anechoic scan robot setup 60Hz-3kHz
date2 = '2021-02-18' #scan 6 & 7 10k-100k AEgir 7 Ran receiving (Motion problems)
date3 = '2021-02-19' #scan 3 maybe best anechoic test 10k-100k lots of clipping on cal signal still
date4 = '2021-02-24' #scan 3,7,8(noise),9(noise)
date = date4

scan = '3'
signal = 'Chirp 10-100kHz 3V'
#desire is the list of all the scans you care to actually look at in this analysis
#desire = [0,4,8] #3 point anechoic x-axis scan
#desire = [0,6,12,18,24] #larger anechoic scan 5 point
desire = [0,12,24]
channels = [1] #recording channels of interest
cal_channel = [1] #recorded calibration of interest (currently should only be len(cal_channel)=1)
legend = ['gen(t)','cal(t)']
#legend = ['Near Wall','Middle','Anechoic']
test = f'Anechoic Test ({scan})'

depth = 0.474 #total depth of the water int he tank
temp = 19
fs = 1e6 #sampling frequency

startzero = 0.5 #leading zeros of the signal
sigL = 2 #signal length
trail = 0.5 #trailing zeros

Acal = (0.6,2.14,depth/2)
Rcal = (0.6,2.06,depth/2)
year = date[0:4]
path = f'W:/uw-measurements-tank/{year}/{date}/{date}_scan'
trec = startzero + sigL + trail #time record length in sec
N = fs*trec
font = 24
###############################################################################
###############################################################################
##############################################################################
##############################################################################
#load generated signal and calibration measurement
gen, cal, ch0, ch1, ch2, ch3 = data.ESAUdata(path+scan, desire, channels, N)

#Load in positions, calculate range distance, plot scan positions desired
#****This will need to be updated for new ESAU version allowing for changing cal
#****and may require need to adjust for plots to not show cal positions

A,R,d = pose.ESAUpose(path+scan, desire, plot=False, Acal = Acal, Rcal = Rcal)

#A = ((0.7,2.1,0.32))
#R = ((0.7,2.02,0.32))
#time delays now allowing for a single measurement w/ single source/receiver pose
if len(A) and len(R) == 3:
    _,_,_,c = tg.timegateTank(A, R, D=depth,T=temp,S=0.03, Coordinate='tank')
else:
    for i in desire:
        _,_,_,c = tg.timegateTank(A[i], R[i], D=depth,T=temp,S=0.03, Coordinate='tank')


##############################################################################
##############################################################################


"""
#Calculate the overall sound pressure level with respect to 1 microPascal as is
#standard for underwater acoustics, rounded to 4 decimal places. 
print('')
print('calculating OASPL re 1e-6 Pa from timewaveform...')
OASPL = np.empty(len(desire))
for i in range(len(desire)):
    OASPL[i] = np.round(10 * np.log10( np.mean( np.square( ch0[:,i] ) ) /1e-6**2),4)
print(f'OASPL{legend} = {OASPL}')
"""

##############################################################################
##############################################################################


print('')
print('plotting waveforms')
print('')
Time = len(cal[:,0])/fs
t = np.arange(0,Time,1/fs)
plt.figure()
for i in range(len(cal_channel)):
    plt.plot(t,gen)
    plt.plot(t,cal[:,i]-np.mean(cal[:,i]))
plt.title(f'{test} {signal} {date}')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (V)')
plt.legend(legend)


##############################################################################
##############################################################################

#need to update MeasGreen to handle various channels.and update the inputs
#to use ch number instead of allsignals. 
print('')
print('')
print('System Response...')
cal1 = np.ndarray.flatten(cal[:,cal_channel]) #obtain only the cal of interest for calculating the system response
hsys,tsys,Hsys,fsys = sys.SystemResponse(cal1,gen,fs,Acal,Rcal,startzero)

#from ESAUResponse import TankResponse
#H_tank,f = TankResponse(ch2,gen,fs,Hsys)





"""
print('')
print('')
print('calculating Measuring Greens function for each point via deconvolution in freq domain...')
H = np.empty((len(gen),len(desire)),dtype = complex)
#for i in range(len(desire)):
for idx,i in enumerate(desire):
    H[:,idx],f = meas.MeasGreen(ch0[:,idx],gen,fs,hsys,A[i],R[i])

#Obtain the single-sided Greens & Freq. Array
Hss = 2*H[0:(int(len(H)/2)),:]
fss = f[0:int(len(f)/2)]
Imp = np.fft.ifft(Hss)

print('')
print('plotting impulse response from deconvolution')
print('')
plt.figure()
for i in range(len(desire)):
    plt.plot(Imp[:,i])
plt.legend(legend)
plt.xlabel('Samples')
plt.ylabel('Amplitude (V)')
plt.title(f'IR {test} {signal}')
#plt.xlim(0,8000)
"""



##############################################################################
##############################################################################
"""
print('')
print('calculating IR & FRF from cross correlation')
print('')
#Impulse Response from the cross correlation for each channel
x0 = np.empty((int(N),len(desire)),dtype = complex)
#Frequency Response (FRF) from fft(time signal)
y0 = np.empty((int(N),len(desire)),dtype = complex)
#single sided FRF
yss = np.empty((int(N/2),len(desire)),dtype = complex)



for i in range(len(desire)):
    xtemp = sci.correlate(ch0[:int(N),i],ch0[:int(N),0],mode='full',method='auto')
    x0[:,i] = xtemp[int(N-1):]
    #to do the fft of the IR properly to obtain an FRF requires 
    #windowing and an inverse filter which has not been done. 
    #so the FRF is calculated using simpy an fft of the time domain. 
    #y0[:,i] = np.fft.fft(x0[:,i])   
    #yss[:,i] = 2*y0[0:(int(len(y0[:,i])/2)),i]
    y0[:,i] = np.fft.fft(ch0[:int(N),i])
    yss[:,i] = 2*y0[0:(int(len(y0[:,i])/2)),i]
    
f = np.fft.fftfreq(len(x0[:,0]),1/fs)
fss = f[:int(len(f)/2)] #single sided freq array
Time = len(x0[:,0])/fs
t = np.arange(0,Time,1/fs)
t = t*1000



#need to EDIT with the below link
#https://www.askpython.com/python/examples/rmse-root-mean-square-error
plt.figure()
for i in range(len(desire)):
    error0 = np.abs(yss[:,i]-yss[:,0])/yss[:,0]*100
    plt.plot(fss,error0)
plt.title('Consistency of Measurements in Same position')
plt.xlabel('Frequency (Hz')
plt.ylabel('% Error')
"""


"""
plt.figure()
for i in range(len(desire)):
    #from byuarglib import psdcalc
    x = ch0[:int(N),i]
    #ns = int(2**np.floor(np.log2(x.size)))
    psd,fpsd,OASPL = byu.psdcalc.psdcalc(x,fs)
    
    plt.plot(fpsd,psd)
plt.title(f'Power Spectral Density of Noise ({legend[i]})')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Square Pressure (Pa^2)')
plt.legend(legend)
"""

"""
print('')
print('calculating IR & FRF from cross correlation')
print('')
#Impulse Response from the cross correlation for each channel
x0 = np.empty((len(gen),len(desire)),dtype = complex)
#Frequency Response (FRF)
y0 = np.empty((len(gen),len(desire)),dtype = complex)
#single sided FRF
yss = np.empty((int(len(gen)/2),len(desire)),dtype = complex)

for i in range(len(desire)):
    xtemp = sci.correlate(ch2[:,i],cal[:,cal_channel[0]],mode='full',method='auto')
    x0[:,i] = xtemp[int(len(gen)-1):]
    y0[:,i] = np.fft.fft(x0[:,i])
    yss[:,i] = 2*y0[0:(int(len(y0[:,i])/2)),i]
f = np.fft.fftfreq(len(x0[:,0]),1/fs)
fss = f[:int(len(f)/2)] #single sided freq array
Time = len(x0[:,0])/fs
t = np.arange(0,Time,1/fs)
t = t*1000
"""


"""
#this is to look at just up to 200kHz since the anechoic panels are rated to 
#only 20k-200k, though the hydrophones are only rated for 100k+ so this includes
#100k-200k bandwidth
idx = np.argwhere(fss==200e3)
idx = idx[0,0]
yssband = np.empty((len(yss[:idx,0]),len(desire)),dtype = complex)
fssband = fss[:idx]
for i in range(len(desire)):
    yssband[:,i] = yss[:idx,i] 
"""

"""
print('')
print('plotting Impulse Response from xcorr..')
print('')
plt.figure()
for i in range(len(desire)):
    plt.plot(t,x0[:,i])
plt.legend(legend)
plt.xlabel('Time (ms)')
plt.ylabel('Amplitude (V)')
plt.title(f'Xcorr IR {test} {signal}')
#plt.xlim(0,8000)


print('')
print('plotting Frequency Response from xcorr IR..')
print('')
plt.figure()
for i in range(len(desire)):
    plt.plot(fss,np.abs(yss[:,i]))
plt.legend(legend)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (V)')
plt.title(f'ssFRF {test} {signal}')
"""



"""
gate_dir = np.argwhere(t<0.34e-3)
gate_dir = max(gate_dir[:,0])
gate_min = np.argwhere(t>0.34e-3)
gate_min = min(gate_min[:,0])
gate_max = np.argwhere(t<0.40e-3)
gate_max = max(gate_max[:,0])
IR_gate_dir = np.real(x0[:gate_dir, 0])
IR_gate_anech = np.real(x0[gate_min:gate_max, 2])
IR_gate_acryl = np.real(x0[gate_min:gate_max, 0])

#OASPL 
OASPL = np.empty((1,len(desire)))
for i in range(len(desire)):
    OASPL[:,i] = byu.OASPLcalc(ch0[:,i])
print(OASPL)
"""

#frequency domain level difference from IR
"""
No worky...
Gxx_dir,flevel,OASPL_dir = byu.autospec(IR_gate_dir,fs,pref=1e-6)
Gxx_anech,flevel,OASPL_anech = byu.autospec(IR_gate_anech,fs,pref=1e-6)
Gxx_acryl,flevel,OASPL_acryl = byu.autospec(IR_gate_acryl,fs,pref=1e-6)
diff_FR_OASPL = OASPL_anech - OASPL_acryl
"""
"""
FR_dir = np.fft.fft(x0[:gate_dir, 0])
FR_anech = np.fft.fft(x0[gate_min:gate_max, 2])
FR_acryl = np.fft.fft(x0[gate_min:gate_max, 0])
level_dir = 10*np.log10(np.sum(np.conj(FR_dir)*FR_dir)/1e-6)
level_anech = 10*np.log10(np.sum(np.conj(FR_anech)*FR_anech)/1e-6)
level_acryl = 10*np.log10(np.sum(np.conj(FR_acryl)*FR_acryl)/1e-6)
diff_level = level_anech - level_acryl
diff_echo = level_anech - level_dir


#time domain from IR OASPL Difference
#OASPL = 10*np.log10(np.sqrt(np.sum(np.conj(FR_dir)*FR_dir))/1e-6)  # I think this needs to be freq domain
oaspl1 = 20*np.log10(np.sqrt(np.mean(IR_gate_dir**2))/(1e-6)**2)
oaspl2 = 20*np.log10(np.sqrt(np.mean(IR_gate_anech**2))/(1e-6)**2)
oaspl3 = 20*np.log10(np.sqrt(np.mean(IR_gate_acryl**2))/(1e-6)**2)
echo_diff = oaspl2 - oaspl3
echo = oaspl2 - oaspl1
"""
"""
##############################################################################
##############################################################################
#determine anechoic effect in frequency domain by finding the decibel difference
#between acrylic and anechoic measure. Use of rolling mean to smooth data
diff_levelband = 20*np.log10(yssband[:,2]/yssband[:,0])
difmeanband = np.mean(diff_levelband)

diff_levelall = 20*np.log10(yss[:,2]/yss[:,0])
difmeanall = np.mean(diff_levelall)
"""






##############################################################################
##############################################################################

"""
#Calculate the overall sound pressure level with respect to 1 microPascal as is
#standard for underwater acoustics, rounded to 4 decimal places. This is
#calculating from the IR instead of the timewaveform as seen above. 
print('')
print('calculating IR OASPL re 1e-6 Pa...')
OASPL = np.empty(len(desire))
for i in range(len(desire)):
    OASPL[i] = np.round(10 * np.log10( np.mean( np.square( x0[:,i] ) ) /1e-6**2),4)
print(f'OASPL{legend} = {OASPL}')
"""

##############################################################################
##############################################################################

"""
plt.figure()
for i in range(len(desire)):
    plt.plot(f,2*np.abs(y0[:,i])**2)
plt.legend(legend)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (V)')
plt.xlim(0,6e5)
plt.title(f'FR {test} {signal}')
#lag, x, line, b = plt.xcorr(allsignals[:,0],allmonitors[:,0])    ***?????***
#lag = sci.correlation_lags(allsignals[:,0].size, allmonitors[:,0].size, mode='same')
#t = lag/fs
#tmax = t[np.argmax(x)]
tmax = len(ch0[:,0])/fs
t = np.arange(0,tmax,1/fs)
idxpeak = np.argmax(x0,axis=0)
tpeak = t[idxpeak]
"""












