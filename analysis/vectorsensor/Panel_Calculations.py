# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:13:48 2021

@author: Corey Dobbs
"""
#%%
import numpy as np
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
import UWIntensityFuncs as uwi
import sys
sys.path.append("D:/uw-acoustics-research/uw-meas-codes/byuarglib/")
import byuarglib as byu
plt.close("all")

fs = 150e3
I_ref = 6.61e-19
sig_cond = .01

filepath_panels = "D:/Box/Underwater Measurements/2021-11-12/"
filepath_nopanels = "D:/Box/Underwater Measurements/2021-11-22/"
p = byu.binfileload(filepath_panels, 'ID', 3, 8)
vx = byu.binfileload(filepath_panels, 'ID', 0, 7)
vy = byu.binfileload(filepath_panels, 'ID', 1, 7)
vz = byu.binfileload(filepath_panels, 'ID', 2, 7)

I3x, fss = uwi.intensityTimeAveraged(p,vx,"x",fs)
I3y, fss = uwi.intensityTimeAveraged(p,vy,"y",fs)
I3z, fss = uwi.intensityTimeAveraged(p,vz,"z",fs)

I3x_act = 10*np.log10(abs(np.real(I3x))/I_ref)
I3x_react = 10*np.log10(abs(np.imag(I3x))/I_ref)
I3y_act = 10*np.log10(abs(np.real(I3y))/I_ref)
I3y_react = 10*np.log10(abs(np.imag(I3y))/I_ref)
I3z_act = 10*np.log10(abs(np.real(I3z))/I_ref)
I3z_react = 10*np.log10(abs(np.imag(I3z))/I_ref)

plt.figure()
plt.plot(fss,I3x_react)
plt.plot(fss,I3x_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('X-direction, with panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')

plt.figure()
plt.plot(fss,I3y_react,fss,I3y_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('Y-direction, with panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')

plt.figure()
plt.plot(fss,I3z_react,fss,I3z_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('Z-direction, with panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')

#NO panels, position 3

p = byu.binfileload(filepath_nopanels, 'ID', 3, 8)
vx = byu.binfileload(filepath_nopanels, 'ID', 0, 7)
vy = byu.binfileload(filepath_nopanels, 'ID', 1, 7)
vz = byu.binfileload(filepath_nopanels, 'ID', 2, 7)

I3x, fss = uwi.intensityTimeAveraged(p,vx,"x",fs)
I3y, fss = uwi.intensityTimeAveraged(p,vy,"y",fs)
I3z, fss = uwi.intensityTimeAveraged(p,vz,"z",fs)

I3x_act = 10*np.log10(abs(np.real(I3x))/I_ref)
I3x_react = 10*np.log10(abs(np.imag(I3x))/I_ref)
I3y_act = 10*np.log10(abs(np.real(I3y))/I_ref)
I3y_react = 10*np.log10(abs(np.imag(I3y))/I_ref)
I3z_act = 10*np.log10(abs(np.real(I3z))/I_ref)
I3z_react = 10*np.log10(abs(np.imag(I3z))/I_ref)

plt.figure()
plt.plot(fss,I3x_react)
plt.plot(fss,I3x_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('X-direction, without panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')

plt.figure()
plt.plot(fss,I3y_react,fss,I3y_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('Y-direction, without panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')

plt.figure()
plt.plot(fss,I3z_react,fss,I3z_act)
plt.xlim([0,5000])
plt.ylim([60,160])
plt.grid()
plt.legend(["Reactive Intensity","Active Intensity"])
plt.title('Z-direction, without panels')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity (dB re 6.61e-19 W/m^2)')