# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 18:56:51 2020

This Function is for loading in ESAU scan positions to plot and use. It searches
for files called scan_positions.txt in the file path. This name convention is 
what is used when using ESAU-Motion-UniversalRobots.  

@author: cvongsaw
"""

def ESAUpose(path,desire = [0],plot = False,Acal = (0.6,2.14,0.3), Rcal = (0.6,2.06,0.3)):
    """
    Parameters
    ----------
    path:   string;
            file path name
    desire: list;
            desired scan positions. Should be input as a list of ordered scans
            this defaults to the first position if not specified otherwise in 
            a list.
    plot:   Boolean {True or False}, optional;
            Default is False which does NOT returns the individual scan plots
            True returns individual scan plots with associated Source and 
            Receiver positions as well as range distance in meters.
    Acal:   Tuple, Optional;
            (x,y,z) position of the AEgir calibration measurement
            Default Acal = (0.6,2.14,0.3)
    Rcal:   Tuple, Optional;
            (x,y,z) position of the Ran calibration measurement 
            Default Rcal = (0.6,2.06,0.3)
    Returns
    -------
    A:      list;
            List of AEgir positions for each individual scan
    R:      list;
            List of Ran positons for each individual scan
    dd:     ndarray;
            range distance between Aegir and Ran positions for the desired 
            scan positions. 

    Notes
    -----
    Author: Cameron Vongsawad
    
    The code runs two options, the first for a scan with only one source and
    receiver position (a single measurement) and the other with more than 
    one (aka an actual "scan")
    
    The code can also check if there is a scan_positions.txt file, if there is 
    none, the code gives back zeros and a warning. 
    
    last modified 2/3/2021
    """

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
    import os.path as check
    
    isFile_gen = check.isfile(path+'/scan_positions.txt')
    if isFile_gen == True:
        print('loading scan positions...')
        ################################################
        #load in scan positions as A(x,y,z) and R(x,y,z)
        ################################################
        pos = np.loadtxt(path+'/scan_positions.txt')
        
        ################################################
        #Load and plot positions of a SINGLE measurement 
        #of single source & receiver positions
        ################################################
        if len(pos) == 24:
            a = (pos[0],pos[1],pos[2])
            r = (pos[12],pos[13],pos[14])
            A = a
            R = r
            ##########################################################    
            #create x,y,and z arrays of scan positions for 3D plotting
            #and highlights positions used for a specific scan
            ##########################################################
            print('organizing scan positions...')    
            xA,yA,zA = np.array([]),np.array([]),np.array([])
            xR,yR,zR = np.array([]),np.array([]),np.array([])
           
            xA = A[0]
            yA = A[1]
            zA = A[2]
            xR = R[0]
            yR = R[1]
            zR = R[2]
            ##################################################################
            #calculate all range distances for the full list of scan positions
            ##################################################################
            d = np.sqrt( (xA-xR)**2 + (yA-yR)**2 + (zA-zR)**2 )
            #Plot Scan Grid
            if plot == True:    
                print('plotting scan position...')
                #plot scan positions  
                from mpl_toolkits import mplot3d 
                scan = plt.figure()
                ax = plt.axes(projection ="3d") 
                ax.scatter3D(xA, yA, zA, color = "green")
                ax.scatter3D(xR, yR, zR,color = "red")
                ax.scatter3D(Acal[0],Acal[1],Acal[2], color = "orange",marker = "^")
                ax.scatter3D(Rcal[0],Rcal[1],Rcal[2], color = "orange",marker = "^")
                ax.scatter3D(A[0],A[1],A[2],color = "blue",marker = 's',linewidths = 10)
                ax.scatter3D(R[0],R[1],R[2],color = "blue",marker = 's',linewidths = 10)
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Y (m)')
                ax.set_zlabel('Z (m)')
                ax.set_xlim(0,1.22)
                ax.set_ylim(0,3.66)
                ax.set_zlim(0,0.91)
                ax.set_title(f'S{A}, R{R}, d={round(d,3)}m')
            #########################################################    
            #Only save range distances for the desired scan positions
            #########################################################    
            dd = d
        
            return A, R, dd
        ####################################################
        #Load and plot positions of a scan of measurements
        #w/ multiple source and multiple receiver positions
        ####################################################
        else: 
            A = []
            a = np.zeros(3)
            R = []
            r = np.zeros(3)
            for i in range(len(pos[:,0])): 
                a = (pos[i,0],pos[i,1],pos[i,2])
                r = (pos[i,12],pos[i,13],pos[i,14])
                A.insert(i,a)
                R.insert(i,r)
            ##########################################################    
            #create x,y,and z arrays of scan positions for 3D plotting
            #and highlights positions used for a specific scan
            ##########################################################
            print('organizing scan positions...')    
            xA,yA,zA = np.array([]),np.array([]),np.array([])
            xR,yR,zR = np.array([]),np.array([]),np.array([])
            for i in range(len(A)):
                xA = np.append(xA,A[i][0])
                yA = np.append(yA,A[i][1])
                zA = np.append(zA,A[i][2])
                xR = np.append(xR,R[i][0])
                yR = np.append(yR,R[i][1])
                zR = np.append(zR,R[i][2])
            ##################################################################
            #calculate all range distances for the full list of scan positions
            ##################################################################
            d = np.empty(len(A))
            for i in range(len(A)):
                d[i] = np.sqrt( (xA[i]-xR[i])**2 + (yA[i]-yR[i])**2 + (zA[i]-zR[i])**2 )
            #Plot Scan Grid
            if plot == True:    
                print('plotting scan positions...')
                for i in desire:
                    #plot scan positions  
                    from mpl_toolkits import mplot3d 
                    scan = plt.figure()
                    ax = plt.axes(projection ="3d") 
                    ax.scatter3D(xA, yA, zA, color = "green")
                    ax.scatter3D(xR, yR, zR,color = "red")
                    ax.scatter3D(Acal[0],Acal[1],Acal[2], color = "orange",marker = "^",linewidths = 4)
                    ax.scatter3D(Rcal[0],Rcal[1],Rcal[2], color = "orange",marker = "^",linewidths = 4)
                    ax.scatter3D(A[i][0],A[i][1],A[i][2],color = "blue",marker = 's',linewidths = 10)
                    ax.scatter3D(R[i][0],R[i][1],R[i][2],color = "blue",marker = 's',linewidths = 10)
                    ax.set_xlabel('X (m)')
                    ax.set_ylabel('Y (m)')
                    ax.set_zlabel('Z (m)')
                    ax.set_xlim(0,1.22)
                    ax.set_ylim(0,3.66)
                    ax.set_zlim(0,0.91)
                    ax.set_title(f'S{A[i]}, R{R[i]}, d={round(d[i],3)}m')
            #########################################################    
            #Only save range distances for the desired scan positions
            #########################################################    
            dd = np.empty(len(desire))
            for idx,i in enumerate(desire): 
                dd[idx] = d[i]
        
            return A, R, dd
    else:
        print('')
        print('Warning: no scan position file found')
        A = [0,0,0]
        R = [0,0,0]
        dd = 0
        return A, R, dd