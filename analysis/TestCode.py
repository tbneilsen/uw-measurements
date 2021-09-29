# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:21:13 2020

TEST CODE

@author: cvong
"""

import numpy as np
#import byuarglib as byu
#import sys as directory
#keep using the natural working directory
#directory.path.insert(0,'C:/Users/cvongsaw/Box/UW Research/Code/uw-measurements')
#add a second working directory
#directory.path.insert(1,'C:/Users/cvongsaw/Box/UW Research/Code/underwater-measurements/analysis/')
from readLogFile import readLogFile
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
date5 = '2021-02-26' #TL measurement [scan 3(5 points) and 4(10 points) good]
date6 = '2021-03-23' #short cal length long signal length
date7 = '2021-05-03' #TL measure (10 points) Scan 1 and 3. Errors with position saving
date8 = '2021-08-03B'#LARGE SCAN for alpha measurement w/out anechoic (great IR not great for T60)
date9 = '2021-08-04' #LARGE SCAN for alpha measurement w/ anechoic (great IR not great for T60)
date10 = '2021-08-06' #Test to improve scans for T60 capability limits
date11 = '2021-08-09'
date12 = 'Need to redo date12 good tests with panels in'
date = date11

scan = '17'
#desire is the list of all the scans you care to actually look at in this analysis
#desire = [0,4,8] #3 point anechoic x-axis scan
#desire = [0,6,12,18,24] #larger anechoic scan 5 point
desire = [0]#,19,56,80]
filename = 'ID000_001log.txt'
calfile = 'cal_000log.txt'
channels = [0,1] #recording channels of interest
cal_channel = [1] #recorded calibration of interest (currently should only be len(cal_channel)=1)
legend = ['0']
#legend = ['Near Wall','Middle','Anechoic']

year = date[0:4]
path = f'W:/uw-measurements-tank/{year}/{date}/{date[0:10]}_scan{scan}/'

freqMin,freqMax,temp,fs,startzero,sigL,trail,totdur,depth,xSource,ySource,zSource,xRec,yRec,zRec = readLogFile(filename,path)
_,_,_,_,_,_,_,_,_,xA_cal,yA_cal,zA_cal,xR_cal,yR_cal,zR_cal = readLogFile(calfile,path)
temp = 19.866
depth = depth
downsample = True

N = fs*sigL #number of samples
signal = f'{sigL}s Chirp {freqMin/1000}kHz-{freqMax/1000}kHz'
test = f'{date} scan{scan} {signal}'


#generated CALIBRATION signal
fscal = fs #sampling frequency
startzerocal = 0.25 #leading zeros of the signal
sigLcal = sigL #signal length
trailcal = 0.25 #trailing zeros
treccal = sigLcal #time record length in sec
Ncal = fscal*treccal
Acal = (xA_cal,yA_cal,zA_cal)#(0.6,2.14,depth/2)
Rcal = (xR_cal,yR_cal,zR_cal)#(0.6,2.06,depth/2)

###############################################################################
###############################################################################
##############################################################################
##############################################################################
#load generated signal and calibration measurement
#gen, calgen, cal, ch0, ch1, ch2, ch3 = data.ESAUdata(path+scan, desire, channels, N, Ncal)

gen, _, cal, ch0, ch1, _, _ = data.ESAUdata(path, desire, channels, N, Ncal)
calgen=gen

#Need to downsample the signal for the T60 code to remove extra noise to IR**2 to at most 2.5x fmax
if downsample == True:
    #resample the data for higher sampling rates???
    factor = 20 #factor by which to reduce signal
    gen = sci.decimate(gen,factor) 
    calgen = sci.decimate(calgen,factor)
    ch0 = sci.decimate(ch0[:,0],factor)
    ch1 = sci.decimate(ch1[:,0],factor)
    fs = fs/factor #"""
    caldec = np.empty([len(ch0),4])
    for i in range(4):
        caldec[:,i] = sci.decimate(cal[:,i],factor)
    cal = caldec
   
#Load in positions, calculate range distance, plot scan positions desired
A,R,dd = pose.ESAUpose(path, desire, plot=False, Acal = Acal, Rcal = Rcal)


#time delays now allowing for a single measurement w/ single source/receiver pose
c = tg.uwsoundspeed(D=depth,T=temp,S=0.03, model='Garrett')


#need to update MeasGreen to handle various channels.and update the inputs
#to use ch number instead of allsignals. 
print('')
print('')
print('System Response...')
cal1 = np.ndarray.flatten(cal[:,cal_channel]) #obtain only the cal of interest for calculating the system response
#cal1 = cal[:,cal_channel[0]] ###!!change!!!####
#cal1 = ch0[:,0]
#ch1 = ch1[:,desire]
#ch1 = ch1[:,0] ##!!change!!!####

tgate,_,tdir,dis = tg.gateValue(A[0],R[0],c) 
tgatec,_,tdirc,dis = tg.gateValue(Acal,Rcal,c) 

#how much to gate the signal in samples
tb4gate = 0.1 #ms before first reflection tgate
Nb4gate = tb4gate/1000 *fs #convert to samples before gating
Ngate = tgate*fs-Nb4gate

tt = np.linspace(0,len(cal1)/fs,len(cal1))      #time array for recorded signal
tt = tt*1000 
"""plt.figure()
plt.plot(tt,ch1)
plt.title(f'Recorded Signal \n {signal} fs={fs/1000}kHz \n {date} scan{scan}')
plt.xlabel('time (ms)')
plt.ylabel('Amplitude')"""

hsys,tsys,Hsys,fsys = response.SysResponse(cal1,calgen,fscal,tgate=tgate,wiener=True,domain='f')
if downsample == False:
    Htank,ftank = response.TankResponse(ch1[:,0],gen,fs,hsys,wiener=True,domain='f')
if downsample == True:
    Htank,ftank = response.TankResponse(ch1,gen,fs,hsys,wiener=True,domain='f')
htank = np.real(np.fft.ifft(Htank))
ttank = np.linspace(0,len(htank)/fs,len(htank))


"""Hss = 2*Hsys[0:(int(len(Hsys)/2))]   #convert to single-sided FRF
fss = fsys[0:int(len(fsys)/2)]/1000    #convert to single-sided from Hz to kHz
Hss_dB = 10*np.log10(np.abs(Hss))   #convert to Levels (dB)

plt.figure()
plt.plot(fss,Hss_dB)
plt.title(f'FRF of {signal} for System')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Amplitude')
plt.grid()

plt.figure()
plt.plot(tsys*1000,hsys)
plt.axvline(tgatec*1000,linestyle='dashed',color='g')
#plt.axvline(Ngate*1000,linestyle='dashdot',color='r') #no idea how to gate this relative to fs
#need to compare this with the fs = 1M and fs = 10M measurements. 
plt.axvline(tdirc*1000,linestyle='dashed',color='r')
plt.xlabel('Time (ms)')
plt.ylabel('Amplitude')
plt.legend(['Chirp IR','Estimated First Reflection','Estimated Direct Signal'])
plt.title(f'Time-Gated Calibration IR of {signal} fs ={fs}Hz')
plt.xlim(0,0.8)
plt.grid()"""





"""___Roll the Shape of the Time-Domain for htank___"""
#The IR is shifted to the end of the array, such that the tail
#spills over to the beginning of the ray. The array must be rolled
#for alignment w/ the actual IR. However, the number of zeros must
#be equal for both leading and trailing zeros. 
roll = True
if roll == True:
    rollt = int(0.5*len(htank)-1)
    shift = rollt/fs
    htank = np.roll(htank,rollt)
else: 
    shift = 0
"""Htss = 2*Htank[0:(int(len(Htank)/2))]   #convert to single-sided FRF
ftss = ftank[0:int(len(ftank)/2)]/1000    #convert to single-sided from Hz to kHz
Htss_dB = 10*np.log10(np.abs(Htss))   #convert to Levels (dB)

plt.figure()
plt.plot(ftss,Htss_dB)
plt.title(f'FRF of {signal} for Tank')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Amplitude')
plt.grid()
#"""

plt.figure()
plt.plot(ttank,htank)
plt.axvline(tgate+shift,linestyle='dashed',color='g')
#plt.axvline(Ngate*1000,linestyle='dashdot',color='r') #no idea how to gate this relative to fs
#need to compare this with the fs = 1M and fs = 10M measurements. 
plt.axvline(tdir+shift,linestyle='dashed',color='r')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.legend(['Chirp IR','Estimated First Reflection','Estimated Direct Signal'])
plt.title(f'Calibrated IR of {signal} fs ={fs}Hz')
#plt.xlim(0,3)
plt.grid()


tbound = tank.T60meas_bounds(htank,fs)
T60 = tank.T60meas(htank,fs,tbound[0],tbound[1],d=depth,c=c,rt='T10',plot=True)

#"""


octData,OctFreq = tank.OctaveFilter(htank,freqMin,freqMax,fs,frac=12)
octTrans = np.transpose(octData)

#propagation absorption estimated for the water characteristics over desired Octave Bands
#prop=0 #when desired to not look into propagation effects. 
prop = tank.alpha_prop(OctFreq,T=temp,S=5,pH=7.2,depth=depth)

    
plt.figure()
legend1 = []
for i in range(len(OctFreq)):
    plt.plot(ttank,octData[i,:])
    legend1.append(np.round(OctFreq[i]))
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('IR Octave Band')
plt.legend(legend1)

plt.figure()
legend = []
for i in range(len(OctFreq)):
    plt.plot(ttank,octTrans[:,i])
    legend.append(np.round(OctFreq[i]))
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('IR Octave Band decimated')
plt.legend(legend1)

octs = 3 #int(len(OctFreq)) #number of octave bands to look at
T60_f = np.empty(int(len(OctFreq[0:octs])))
a_wall_f = np.empty(int(len(OctFreq[0:octs])))
for i in range(len(OctFreq[0:octs])):
    tbound = tank.T60meas_bounds(octTrans[:,i],fs)
    T60_f[i] = tank.T60meas(octTrans[:,i],fs,tbound[0],tbound[1],d=depth,c=c,rt='T10',plot=True)
    a_wall_f[i] = tank.alpha_wall(T60_f[i],d=depth,c=c,acc=False,anech=False,alpha_p=0)#prop[i])
#"""

#overall absorption coefficient over entire bandwidth
a_wall_gen = tank.alpha_wall(T60,d=depth,c=c,acc=False,anech=False,alpha_p=0)#prop)

T60,sig,fschroeder = tank.T60est(depth,c=c,)

"""
ht = np.abs(np.fft.ifft(Htank))
T60 = tank.T60meas(ht,fs,t0=0,t1=int(len(ht)/fs),d=depth,c=c,rt='T60',plot=True)
alpha_s = tank.alpha_wall(T60,depth,c,acc=False,anech=False,alpha_p=0)
"""
"""
Hss = 2*Hsys[0:(int(len(Hsys)/2))]   #convert to single-sided FRF
fss = fsys[0:int(len(fsys)/2)]/1000    #convert to single-sided from Hz to kHz
Hss_dB = 10*np.log10(Hss) 

Hsstank = 2*Htank[0:(int(len(Htank)/2))]   #convert to single-sided FRF
fsstank = ftank[0:int(len(ftank)/2)]/1000    #convert to single-sided from Hz to kHz
Hsstank_dB = 10*np.log10(Hsstank/1e-6) 
"""
"""
ir = response.IR(ch1,gen,fs,wiener=False,domain='f')
ir1 = 10*np.log10(ir**2)
tt = np.linspace(0,len(ir1)/fs,len(ir1))               #time array for ht
tt = tt*1000 
frf = np.fft.fft(ir)
frf = 10*np.log10(frf/1e-6)
frf = 2*frf[0:(int(len(frf)/2))]
f = np.fft.fftfreq(len(ir),d=1/fs)
f = f[0:int(len(f)/2)]/1000
plt.figure()
plt.plot(tt,ir)
plt.title('straight IR')
plt.xlabel('time (ms)')
plt.ylabel('')
plt.grid()
plt.figure()
plt.plot(f,frf)
plt.title(f'Straight FRF')
plt.xlabel('Frequency (kHz)')
plt.ylabel(r'Level (dB re 1 $\mu Pa$)')
plt.grid()
"""




"""
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
"""







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














