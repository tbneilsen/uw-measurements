"""
Created on Thurs Jun 2 15:11:18 2022

@author: Alexandra Hopps

"""
import os
import matplotlib
import math as m
from seabird.cnv import fCNV 
from matplotlib import pyplot as plt
import numpy as np


def plot_CTD(profile_path, save_path, fig_type, filename_adj = '', scale = 5): 
    '''
    This function plots the sound speed profiles of CTD data in .cnv files. 
    The function is designed to iterate through all the data in a folder and save a plot for each CTD file. 
    Additional changes can be made within the function to edit the plotting parameters. 
    
    Inputs
    ------
    profile_path     : string; 
                       This is the path to a folder with the CTD data and should end with a slash.
    save_path        : string; 
                       This is the path that all the figures will be saved to should end with a slash.
    fig_type         : string; 
                       The type of figure you want to save the figure as (.jpg or .png)
    filename_adj     : string, optional; 
                       Any adjustments you want to make to the filename to make it more unique. 
                       Defaults to an empty string. 
                       Example of using filename_adj: filename_adj = 'DayOne' and now the name of the file is <CTD#>SoundSpeedProfileDayOne.jpg instead of <CTD#>SoundSpeedProfile.jpg 
                                                       
    Outputs
    -------
    Returns a Sound Speed Profile for each file in the profile_path folder. 

    NOTE to see the keys from the CTD .cnv file, run the following code:
        print("Header: %s" % profile.attributes.keys())
        print("Data: %s" % profile.keys()) 
    You can use this to obtain more information about the CTD and plot different parameters (e.g. Depth vs Temperature).

    '''
   
   # Update the parameters of the plot
    params = {'legend.fontsize': 25,
              'figure.figsize': (20, 15),
             'axes.labelsize':32,
             'axes.titlesize':35,
             'axes.titleweight':'bold',
             'xtick.labelsize':26,
             'ytick.labelsize':26,
             'lines.linewidth':1}
    matplotlib.rcParams.update(params)

    files = os.listdir(profile_path) # Makes an array of all files in the CTDD data folder

    for file in files:
        save_name = file[5:-4] + 'SoundSpeedProfile.jpg' 
        profile = fCNV(profile_path + file)  

        plt.figure() 
        plt.plot(profile['soundspeed'], profile['DEPTH'])
        plt.title(file[5:-4]+'Sound Speed Profile')
        plt.ylabel('Depth (m)')
        plt.xlabel('Sound Speed (m/s)')
        min_x = m.floor((min(profile['soundspeed']) - 0.05) * scale)/scale # Makes the axis of your plot fit the data set
        max_x = m.ceil((max(profile['soundspeed']) + 0.05) * scale)/scale 

        plt.xlim(min_x,max_x)
        plt.grid()
        plt.savefig(save_path+save_name)
    
# NOTE This code can be easily edited to put all the sound speed profiles on one figure by taking the plt.figue() and plt.savefig() outside of the loop.

profile_path = '/home/byu.local/hoppsa/underwater/data/ocean_ssps/sbcex2022-leg2/Data/Test Function Data/'
save_path = '/home/byu.local/hoppsa/underwater/data/ocean_ssps/sbcex2022-leg2/Temperature Plots/Test Temperature Plots/'
fig_type = '.jpg'
filename_adj = ''

plot_CTD(profile_path, save_path, fig_type, filename_adj)