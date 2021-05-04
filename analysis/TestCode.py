# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:21:13 2020

TEST CODE

@author: cvong
"""

import numpy as np
import byuarglib as byu
import sys as directory
#keep using the natural working directory
directory.path.insert(0,'C:/Users/cvongsaw/Box/UW Research/Code/uw-measurements')
#add a second working directory
directory.path.insert(1,'C:/Users/cvongsaw/Box/UW Research/Code/underwater-measurements/analysis/')
from readLogFile import readLogFile
import SystemFResponse as sys
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
import ESAUResponse as response
import TankCharacterization as tank
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
date5 = '2021-02-26' # TL measurement [scan 3(5 points) and 4(10 points) good]
date6 = '2021-03-23' #short cal length long signal length
date7 = '2021-05-03' #TL measure (10 points) Scan 1 and 3. Errors with position saving
date = date7

scan = '3'
signal = '0.5s Chirp 10Hz-10kHz'
#desire is the list of all the scans you care to actually look at in this analysis
#desire = [0,4,8] #3 point anechoic x-axis scan
#desire = [0,6,12,18,24] #larger anechoic scan 5 point
desire = [0]#,1,2,3,4,5,6,7,8,9]
channels = [1] #recording channels of interest
cal_channel = [1] #recorded calibration of interest (currently should only be len(cal_channel)=1)
legend = ['0']#,'1','2','3','4','5','6','7','8','9']
#legend = ['Near Wall','Middle','Anechoic']
test = f'{date3} scan{scan} {signal}'

year = date[0:4]
path = f'W:/uw-measurements-tank/{year}/{date}/{date}_scan'
path2 = f'W:/uw-measurements-tank/{year}/{date}/{date}_scan{scan}/'
filename = 'ID000_001log.txt'

freqMin,freqMax,temp,fs,sigL,depth,xSource,ySource,zSource,xRec,yRec,zRec = readLogFile(filename,path2)
startzero = 0.5 #leading zeros of the signal
trail = 2.49996 #trailing zeros
trec = startzero + sigL + trail #time record length in sec
N = fs*trec

#generated CALIBRATION signal
fscal = fs #sampling frequency
startzerocal = 0.5 #leading zeros of the signal
sigLcal = sigL #signal length
trailcal = 0.5 #trailing zeros
treccal = startzerocal + sigLcal + trailcal #time record length in sec
Ncal = fscal*treccal
Acal = (0.6,2.14,depth/2)
Rcal = (0.6,2.06,depth/2)

###############################################################################
###############################################################################
##############################################################################
##############################################################################
#load generated signal and calibration measurement
#gen, calgen, cal, ch0, ch1, ch2, ch3 = data.ESAUdata(path+scan, desire, channels, N, Ncal)

gen, _, cal, _, ch1, _, _ = data.ESAUdata(path+scan, desire, channels, N, Ncal)
#calgen = gen


#Load in positions, calculate range distance, plot scan positions desired
A,R,dd = pose.ESAUpose(path+scan, desire, plot=False, Acal = Acal, Rcal = Rcal)

#A = ((0.7,2.1,0.32))
#R = ((0.7,2.02,0.32))
#time delays now allowing for a single measurement w/ single source/receiver pose
c = tg.uwsoundspeed(D=depth,T=temp,S=0.03, model='Garrett')




##############################################################################
"""Plot Waveforms"""
##############################################################################

"""
print('')
print('plotting waveforms')
print('')
Time = len(cal[:,0])/fs
t = np.arange(0,Time,1/fs)
plt.figure()
for i in range(len(desire)):
    #plt.plot(t,gen)
    plt.plot(t,ch1[:,i]-np.mean(ch1[:,i]))  # 0 mean to remove DC offset
plt.title(f'{test} {signal} {date}')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (V)')
plt.legend(legend)
"""
##############################################################################
##############################################################################

#need to update MeasGreen to handle various channels.and update the inputs
#to use ch number instead of allsignals. 
print('')
print('')
print('System Response...')
#cal1 = np.ndarray.flatten(cal[:,cal_channel]) #obtain only the cal of interest for calculating the system response
cal1 = cal[:,cal_channel[0]] ###!!change!!!####
#ch1 = ch1[:,desire]
ch1 = ch1[:,0] ##!!change!!!####

tgate,_,dis = tg.gateValue(Acal,Rcal,c) 
tdir = dis/c
#how much to gate the signal in samples
tb4gate = 0.1 #ms before first reflection tgate
Nb4gate = tb4gate/1000 *fs #convert to samples before gating
Ngate = tgate*fs-Nb4gate

hsys,tsys,Hsys,fsys = sys.SystemResponse(cal1,gen,fscal,Acal,Rcal,tgate,bandpass = False,wiener=True)
#hsys1,tsys1,Hsys1,fsys1 = sys.SystemResponse(ch1,gen,fscal,Acal,Rcal,tgate,bandpass = True,wiener=True)

Htank,ftank = response.TankResponse(ch1,gen,fs,hsys,bandpass=True,wiener=True)


ht = np.abs(np.fft.ifft(Htank))
T60,alpha_S = tank.T60meas(ht,fs,depth,c,rt='T60',plot=True,acc=False,alpha_p=0)


Hss = 2*Hsys[0:(int(len(Hsys)/2))]   #convert to single-sided FRF
fss = fsys[0:int(len(fsys)/2)]/1000    #convert to single-sided from Hz to kHz
Hss_dB = 10*np.log10(Hss) 

#Hss1 = 2*Hsys1[0:(int(len(Hsys1)/2))]   #convert to single-sided FRF
#fss1 = fsys1[0:int(len(fsys1)/2)]/1000    #convert to single-sided from Hz to kHz
#Hss1_dB = 10*np.log10(Hss1/1e-6) 

Hsstank = 2*Htank[0:(int(len(Htank)/2))]   #convert to single-sided FRF
fsstank = ftank[0:int(len(ftank)/2)]/1000    #convert to single-sided from Hz to kHz
Hsstank_dB = 10*np.log10(Hsstank/1e-6) 

plt.figure()
#plt.plot(tsys,hsys)
plt.plot(tsys,hsys)
plt.axvline(tgate*1000,linestyle='dashed',color='g')
#plt.axvline(Ngate*1000,linestyle='dashdot',color='r') #no idea how to gate this relative to fs
#need to compare this with the fs = 1M and fs = 10M measurements. 
plt.axvline(tdir*1000,linestyle='dashed',color='r')
plt.xlabel('Time (ms)')
plt.ylabel('Amplitude')
plt.legend(['Chirp IR','Estimated First Reflection','Estimated Direct Signal'])
plt.title(f'Time-Gated Calibration IR of {signal} fs ={fs}Hz')
plt.xlim(0,3)
plt.grid()

plt.figure()
plt.plot(fss,Hss_dB)
plt.title(f'Time-Gated Calibration FRF of {signal} Swept-Sine')
plt.xlabel('Frequency (kHz)')
plt.ylabel(r'Level (dB re 1 $\mu Pa$)')
#plt.xlim(0,300)
#plt.ylim(35,65)
plt.grid()

plt.figure()
plt.plot(fsstank,Hsstank_dB)
plt.title(f'Calibrated Tank Transfer Function of {signal} Swept-Sine')
plt.xlabel('Frequency (kHz)')
plt.ylabel(r'Level (dB re 1 $\mu Pa$)')
#plt.xlim(0,300)
#plt.ylim(35,65)
plt.grid()


#explore if this system response in reverse can give the recorded value
#invest = sci.correlate(hsys1,np.flip(gen),mode='full',method='auto')
#plt.figure()
#plt.plot(ch1)
#plt.plot(invest)
#plt.title('recorded measure and calculated measure (from time domain)')


"""
#need to EDIT with the below link
#https://www.askpython.com/python/examples/rmse-root-mean-square-error

"""







##############################################################################
##############################################################################
#OASPL w/ and w/out panels
##############################################################################
##############################################################################

##############################################################################
"""OASPL"""
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

