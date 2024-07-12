# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:05:25 2024

@author: Carlo Vignali

@def: This is preprocessing code for fiber photometry data

@source: Modified from FibPhoEpocAveraging.py written by David Root (in sourceCode folder)
"""
#%% Part 0.0 - import packages

# special call for jupyter notebook only
#get_ipython().run_line_magic('matplotlib','inline')

import os
import numpy as np # fundamental package for scientific computing, handles arrays and maths
from sklearn.metrics import auc
import matplotlib.pyplot as plt  # standard Python plotting library
import scipy.stats as stats
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import linregress
from scipy.optimize import curve_fit, minimize

import tdt # tdt library

#

os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
import localFunctions as lf

#%% Part 0.1 - logical chain

if not 'userInputs' in locals():
# !! SET LOGIC CHAIN !!
    print('These will be logical decisions, so type either 1 or 0')
    userInputs = {'doBxOnly':[], 'pullMPC':[], 'pullSynapse':[], 'errorPullMPC':[]}
    userInputs['doRxOnly'] = bool(int(input('Do you want to analyze ONLY the recording data?: ')))
    userInputs['oneMouse'] = bool(int(input('Do you want to analyze ONE mouse [alt. is group analysis]?: ')))
    if not userInputs['doBxOnly']:
        userInputs['pullSynapse'] = bool(int(input('Do you have fiber photometry data associated with this MPC file?: ')))
        userInputs['pullMPC'] = bool(int(input('Do you need to preprocess the MPC file(s)?: ')))
    else:
        userInputs['pullMPC'] = True    
    
    # !! SET PATH !!
    inputPath = ('S:\\Private\\Data\\Vignali cv105\\data\\MARS\\')
    outputPath = ('S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\')
    
    if userInputs['oneMouse']:
        # !! SET MOUSE !!
        #determine path to datafile - experiment file + animal number + date
        currFolders = os.listdir(inputPath)[1:] #from path, list folders
        print(currFolders)
        # INPUT THE MOUSE ID OF INTEREST AND NOT THE PYTHON INDEX FOR THAT MOUSE
        thisMouse = currFolders[int(input('Which mouse from currFolders: '))-1] #then user provides input to choose the animal of interest
    

#%% Part 0.2 - format data folder and set root path

if userInputs['oneMouse']:
    # !! MOVE FILES !!
    lf.moveFiles(thisMouse,inputPath,suffix='-')
    if date == []:
       date = input('Date not found: input experiment date [YYMMDD]: ') 
       
    # !! SET PATH !!
    dataPath = inputPath + thisMouse + '\\' + date + '\\'
    bxPath = outputPath + thisMouse + '\\bx\\'
    rxPath = outputPath + thisMouse + '\\' + date + '\\fiber\\'
    os.makedirs(os.path.dirname(rxPath), exist_ok=True)

#%% Part 1.0 - run MedPC analysis script

    if not userInputs['doRxOnly']:
        if userInputs['pullMPC']:
            import Parent
            Parent.py
            
        else:
            import pickle
            procDataSheet = pickle.load(open(os.path.join(bxPath, 'suppDataFile.pkl'), 'rb'))
     
#%% Part 1.1 - load data

            # !! LOAD DATATANK !!
            tankPath = dataPath + next(file for file in os.listdir(dataPath) if '{m}-'.format(m=thisMouse) in file) + '\\'
            dataTank = tdt.read_block(tankPath, export='interlaced', outdir=tankPath, prefix='iso') #read length of data folder (tankPath)
                # note: the export and outdir flags allow the read_block function to pull out individual data streams for each sensor, and saves files in tankPath as F32 files (these are 32 bit floating point files - commonly used to represent sound waves)
                # ! see print(tdt.read_block.__doc__) for more details !
    

#%% Part 1.2 - import photometry stream

            # !! PULL CHANNEL DATA !!
            rxSite = input('What hemisphere is recording data from [left / right / both]: ')
            if rxSite == 'left':
                data405_LH = np.fromfile(tankPath + 'iso__405A.f32', dtype=np.float32)
                data465_LH = np.fromfile(tankPath + 'iso__465A.f32', dtype=np.float32)
                data560_LH = np.fromfile(tankPath + 'iso__560B.f32', dtype=np.float32)
                data405_RH = np.array([])
                data465_RH = np.array([])
                data560_RH = np.array([])
            elif rxSite == 'right':
                data405_LH = np.array([])
                data465_LH = np.array([])
                data560_LH = np.array([])
                data405_RH = np.fromfile(tankPath + 'iso__405C.f32', dtype=np.float32)
                data465_RH = np.fromfile(tankPath + 'iso__465C.f32', dtype=np.float32)
                data560_RH = np.fromfile(tankPath + 'iso__560D.f32', dtype=np.float32)
            elif rxSite == 'both':
                data405_LH = np.fromfile(tankPath + 'iso__405A.f32', dtype=np.float32)
                data465_LH = np.fromfile(tankPath + 'iso__465A.f32', dtype=np.float32)
                data560_LH = np.fromfile(tankPath + 'iso__560B.f32', dtype=np.float32)
                data405_RH = np.fromfile(tankPath + 'iso__405C.f32', dtype=np.float32)
                data465_RH = np.fromfile(tankPath + 'iso__465C.f32', dtype=np.float32)
                data560_RH = np.fromfile(tankPath + 'iso__560D.f32', dtype=np.float32)
            
            dataName_indivFP = ['data405_LH','data465_LH','data560_LH','data405_RH','data465_RH','data560_RH']
            colors_indivPlot = ['#0000FF','#008000','#FF0000','#ADD8E6','#90EE90','#FFC0CB']
            for f in range(len(dataName_indivFP)):
                if eval(dataName_indivFP[f]).any(): 
                    plt.figure(figsize=(10, 4))
                    plt.plot(eval(dataName_indivFP[f])[500:], color=colors_indivPlot[f])
                    plt.title('Time Series for {}'.format(dataName_indivFP[f]))
                    plt.xlabel('approx. time (ms)')
                    plt.ylabel('raw fluorescence')
                    plt.show()    
    
#%% Part 2.0 - preprocessing section [adapted from 'photometry data preprocessing.ipynb' in sourceCode folder]
# Steps are as follows:
    # 1. Lowpass filtering to reduce noise.
    # 2. Correction for photobleaching, i.e. the slow decreace in the fluorescence signal over time. Two different methods are shown (i) subtraction of a double exponential fit (ii) highpass filtering with a very low cutoff frequency.
    # 3. Movement correction by subtracting a linear fit of the movement control channel.
    # 4. Conversion of the signal to dF/F.
    
    #set default plot properties
            plt.rcParams['figure.figsize'] = [14, 12] # Make default figure size larger.
            plt.rcParams['axes.xmargin'] = 0          # Make default margin on x axis zero.
            plt.rcParams['axes.labelsize'] = 12     #Set default axes label size 
            plt.rcParams['axes.titlesize']=15
            plt.rcParams['axes.titleweight']='heavy'
            plt.rcParams['ytick.labelsize']= 10
            plt.rcParams['xtick.labelsize']= 10
            plt.rcParams['legend.fontsize']=12
            plt.rcParams['legend.markerscale']=2
            
#%% Part 2.1 - 

isos_raw = data405_RH
gACh_raw = data465_RH
rdLight_raw = data560_RH
time_seconds = epoc_ticker
sampling_rate = 1017
plotPoints = np.linspace(0,len(time_seconds),len(rdLight_raw))

##plot raw

#RdLight1 raw
fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
plot1=ax1.plot(plotPoints, rdLight_raw, 'r', label='rdLight') #plot rdLight on left y-axis
ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
plot2=ax2.plot(plotPoints, isos_raw, 'b', label='isos') # plot isos on right y-axis

# Plot rewards times as ticks.
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), 
#                        label='Reward Cue', color='w', marker="|", mec='k')

ax1.set_ylim(np.min(rdLight_raw[500:]), np.max(rdLight_raw[500:]))
ax2.set_ylim(np.min(isos_raw[500:]), np.max(isos_raw[500:]))
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='r')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Raw signals')

lines = plot1 + plot2 #+reward_ticks #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.93)) #add legend

#grabACh3.1+ raw
fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
plot1=ax1.plot(plotPoints, gACh_raw, 'g', label='grabACh') # plot isos on right y-axis
ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
plot2=ax2.plot(plotPoints, isos_raw, 'b', label='isos') # plot isos on right y-axis

# Plot rewards times as ticks.
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), 
#                        label='Reward Cue', color='w', marker="|", mec='k')

ax1.set_ylim(np.min(gACh_raw[500:]), np.max(gACh_raw[500:]))
ax2.set_ylim(np.min(isos_raw[500:]), np.max(isos_raw[500:]))
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('grabACh Signal (V)', color='g')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Raw signals')

lines = plot1 + plot2 #+reward_ticks #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.93)) #add legend

#%% Part 2.2 - denoise with lowpass filter

# Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
b,a = butter(2, 10, btype='low', fs=sampling_rate)
rdLight_denoised = filtfilt(b,a, rdLight_raw)
gACh_denoised = filtfilt(b,a, gACh_raw)
isos_denoised = filtfilt(b,a, isos_raw)

#RdLight1 denoised
fig,ax1=plt.subplots()
plot1=ax1.plot(plotPoints, rdLight_denoised, 'r', label='rdLight denoised')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_denoised, 'b', label='isos denoised')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k', ms=10)

ax1.set_ylim(np.min(rdLight_denoised[500:]), np.max(rdLight_denoised[500:]))
ax2.set_ylim(np.min(isos_denoised[500:]), np.max(isos_denoised[500:]))
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='r')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals')

lines = plot1+plot2 #+reward_ticks #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.92)) #add legend

#grabACh denoised
fig,ax1=plt.subplots()
plot1=ax1.plot(plotPoints, gACh_denoised, 'g', label='grabACh denoised')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_denoised, 'b', label='isos denoised')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k', ms=10)

ax1.set_ylim(np.min(gACh_denoised[500:]), np.max(gACh_denoised[500:]))
ax2.set_ylim(np.min(isos_denoised[500:]), np.max(isos_denoised[500:]))
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('grabACh Signal (V)', color='g')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals')

lines = plot1+plot2 # +reward_ticks #line handle for legend
labels = [l.get_label() for l in lines]  #get legend labels
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.92)) #add legend


#%% Part 2.3 - plot smaller window of post-denoise data

#RdLight1 denoised (small window)

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, rdLight_raw, color='r', alpha=0.3, label='rdLight raw')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_raw, color='b', alpha=0.3, label='isos raw') 
plot3=ax1.plot(plotPoints, rdLight_denoised, color='r', label='rdLight denoised') 
plot4=ax2.plot(plotPoints, isos_denoised, color='b', label='isos denoised') 
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.59), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='r')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals')

lines = plot1+plot2 + plot3 + plot4 #+ reward_ticks
labels = [l.get_label() for l in lines]
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
ax1.set_xlim(2000, 2060) # 60 sec window
ax1.set_ylim(np.min(rdLight_raw[500:]), np.max(rdLight_raw[500:]))
ax2.set_ylim(np.min(isos_raw[500:]), np.max(isos_raw[500:]))

#grabACh denoised (small window)

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, gACh_raw, color='g', alpha=0.3, label='grabAch raw')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_raw, color='b', alpha=0.3, label='isos raw') 
plot3=ax1.plot(plotPoints, gACh_denoised, color='g', label='grabAch denoised') 
plot4=ax2.plot(plotPoints, isos_denoised, color='b', label='isos denoised') 
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.59), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('grabAch Signal (V)', color='g')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals')

lines = plot1+plot2 + plot3 + plot4 #+ reward_ticks
labels = [l.get_label() for l in lines]
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
ax1.set_xlim(1000, 1060) # 60 sec window
ax1.set_ylim(86,96)#np.min(gACh_raw[500:]), np.max(gACh_raw[500:]))
ax2.set_ylim(60,70)#np.min(isos_raw[500:]), np.max(isos_raw[500:]))

#denoised comp

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, gACh_raw, color='g', alpha=0.3, label='grabAch raw')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, rdLight_raw, color='r', alpha=0.3, label='rdLight raw') 
plot3=ax1.plot(plotPoints, gACh_denoised, color='g', label='grabAch denoised') 
plot4=ax2.plot(plotPoints, rdLight_denoised, color='r', label='rdLight denoised') 
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.59), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('grabAch Signal (V)', color='g')
ax2.set_ylabel('rdLight Signal (V)', color='r')
ax1.set_title('Denoised signals')

lines = plot1+plot2 + plot3 + plot4 #+ reward_ticks
labels = [l.get_label() for l in lines]
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
ax1.set_xlim(1000, 1060) # 60 sec window
ax1.set_ylim(86,96)#np.min(gACh_raw[500:]), np.max(gACh_raw[500:]))
ax2.set_ylim(14,24)#np.min(rdLight_raw[500:]), np.max(rdLight_raw[500:]))

#%% Part 2.4 - photobleach correction option #1 > double exponential

# The double exponential curve we are going to fit.
def double_exponential(t, const, amp_fast, amp_slow, tau_slow, tau_multiplier):
    '''Compute a double exponential function with constant offset.
    Parameters:
    t       : Time vector in seconds.
    const   : Amplitude of the constant offset. 
    amp_fast: Amplitude of the fast component.  
    amp_slow: Amplitude of the slow component.  
    tau_slow: Time constant of slow component in seconds.
    tau_multiplier: Time constant of fast component relative to slow. 
    '''
    tau_fast = tau_slow*tau_multiplier
    return const+amp_slow*np.exp(-t/tau_slow)+amp_fast*np.exp(-t/tau_fast)

## !! correction for rdLight !!

# Fit curve to rdLight signal.
max_sig = np.max(rdLight_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
rdLight_parms, parm_cov = curve_fit(double_exponential, plotPoints, rdLight_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
rdLight_expfit = double_exponential(plotPoints, *rdLight_parms)

# Fit curve to isos signal.
max_sig = np.max(isos_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
isos_parms, parm_cov = curve_fit(double_exponential, plotPoints, isos_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
isos_expfit = double_exponential(plotPoints, *isos_parms)

#plot fits over denoised data
fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, rdLight_denoised, 'r', label='rdLight')
plot3=ax1.plot(plotPoints, rdLight_expfit, color='k', linewidth=1.5, label='Exponential fit') 
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_denoised, color='b', label='isos') 
plot4=ax2.plot(plotPoints, isos_expfit,color='k', linewidth=1.5) 


ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='r')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals with double exponential fits')

lines = plot1 + plot2 + plot3
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
ax1.set_ylim(np.min(rdLight_raw[500:]), np.max(rdLight_raw[500:]))
ax2.set_ylim(np.min(isos_raw[500:]), np.max(isos_raw[500:]))

##subtract exponentials fit to data

rdLight_detrended = rdLight_denoised - rdLight_expfit
isos_detrended = isos_denoised - isos_expfit

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, rdLight_detrended, 'r', label='rdLight')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_detrended, color='b', label='isos') 

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='r')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Bleaching Correction by Double Exponential Fit')

lines = plot1+plot2 
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
#ax1.set_xlim(1000, 1060) # 60 sec window
ax1.set_ylim(-3+np.min(rdLight_detrended[500:]), np.max(rdLight_detrended[500:]))
ax2.set_ylim(np.min(isos_detrended[500:]), 3+np.max(isos_detrended[500:]))


## !! correction for gACh !!

# Fit curve to gACh signal.
max_sig = np.max(gACh_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
gACh_parms, parm_cov = curve_fit(double_exponential, plotPoints, gACh_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
gACh_expfit = double_exponential(plotPoints, *gACh_parms)

# Fit curve to isos signal.
max_sig = np.max(isos_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
isos_parms, parm_cov = curve_fit(double_exponential, plotPoints, isos_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
isos_expfit = double_exponential(plotPoints, *isos_parms)

#plot fits over denoised data
fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, gACh_denoised, 'g', label='gACh')
plot3=ax1.plot(plotPoints, gACh_expfit, color='k', linewidth=1.5, label='Exponential fit') 
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_denoised, color='b', label='isos') 
plot4=ax2.plot(plotPoints, isos_expfit,color='k', linewidth=1.5) 


ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('gACh Signal (V)', color='g')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Denoised signals with double exponential fits')

lines = plot1 + plot2 + plot3
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
ax1.set_ylim(np.min(gACh_raw[500:]), np.max(gACh_raw[500:]))
ax2.set_ylim(np.min(isos_raw[500:]), np.max(isos_raw[500:]))

##subtract exponentials fit to data

gACh_detrended = gACh_denoised - gACh_expfit
isos_detrended = isos_denoised - isos_expfit

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, gACh_detrended, 'g', label='gACh')
ax2=plt.twinx()
plot2=ax2.plot(plotPoints, isos_detrended, color='b', label='isos') 

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('gACh Signal (V)', color='g')
ax2.set_ylabel('isos Signal (V)', color='b')
ax1.set_title('Bleaching Correction by Double Exponential Fit')

lines = plot1+plot2 
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right'); 
ax1.set_xlim(1000, 1060) # 60 sec window
ax1.set_ylim(np.min(gACh_detrended[500:]), 3+np.max(gACh_detrended[500:]))
ax2.set_ylim(-3+np.min(isos_detrended[500:]), np.max(isos_detrended[500:]))

#%% Part 2.5 - motion correction

## !! for RdLight1 signal !!

slope, intercept, r_value, p_value, std_err = linregress(x=isos_detrended, y=rdLight_detrended)

plt.scatter(isos_detrended[::5], rdLight_detrended[::5],alpha=0.1, marker='.')
x = np.array(plt.xlim())
plt.plot(x, intercept+slope*x)
plt.xlabel('isos')
plt.ylabel('rdLight')
plt.title('isos - rdLight correlation.')

print('Slope    : {:.3f}'.format(slope))
print('R-squared: {:.3f}'.format(r_value**2))

##plot

rdLight_est_motion = intercept + slope * isos_detrended
rdLight_corrected = rdLight_detrended - rdLight_est_motion

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, rdLight_detrended, 'k' , label='rdLight - pre motion correction', alpha=0.5)
plot3=ax1.plot(plotPoints, rdLight_corrected, 'r', label='rdLight - motion corrected', alpha=0.5)
ax2=plt.twinx()
plot4=ax2.plot(plotPoints, rdLight_est_motion - 0.05, 'y', label='estimated motion')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 0.08), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('rdLight Signal (V)', color='g')
ax1.set_title('Motion Correction')

lines = plot1+plot3+plot4 #+ reward_ticks
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

ax1.set_xlim(1000, 1020)  # 60 sec window
ax1.set_ylim(np.min(rdLight_corrected[500:]), np.max(rdLight_corrected[500:]))
ax2.set_ylim(1+np.min(rdLight_corrected[500:]), 1+np.max(rdLight_corrected[500:]))


## !! for grabACh signal !!

slope, intercept, r_value, p_value, std_err = linregress(x=isos_detrended, y=gACh_detrended)

plt.scatter(isos_detrended[::5], gACh_detrended[::5],alpha=0.1, marker='.')
x = np.array(plt.xlim())
plt.plot(x, intercept+slope*x)
plt.xlabel('isos')
plt.ylabel('gACh')
plt.title('isos - gACh correlation.')

print('Slope    : {:.3f}'.format(slope))
print('R-squared: {:.3f}'.format(r_value**2))

##plot

gACh_est_motion = intercept + slope * isos_detrended
gACh_corrected = gACh_detrended - gACh_est_motion

fig,ax1=plt.subplots()  
plot1=ax1.plot(plotPoints, gACh_detrended, 'k' , label='gACh - pre motion correction', alpha=0.5)
plot3=ax1.plot(plotPoints, gACh_corrected, 'g', label='gACh - motion corrected', alpha=0.5)
ax2=plt.twinx()
plot4=ax2.plot(plotPoints, gACh_est_motion - 0.05, 'y', label='estimated motion')
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 0.08), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('gACh Signal (V)', color='g')
ax1.set_title('Motion Correction')

lines = plot1+plot3+plot4 #+ reward_ticks
labels = [l.get_label() for l in lines]  
legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

ax1.set_xlim(1000, 1060)  # 60 sec window
ax1.set_ylim(-2+np.min(gACh_corrected[500:]),4)# np.max(gACh_corrected[500:]))
ax2.set_ylim(1+np.min(gACh_corrected[500:]), 1+np.max(gACh_corrected[500:]))

#%%normalize option #1

rdLight_dF_F = 100*rdLight_corrected/rdLight_expfit
gACh_dF_F = 100*gACh_corrected/gACh_expfit

fig,axs=plt.subplots(2,1)  
plot1=axs[0,0].plot(plotPoints, rdLight_dF_F, 'r', label='rdLight dF/F')
#reward_ticks = axs1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 6), label='Reward Cue',color='w', marker="v", mfc='k', mec='k', ms=8)

axs[0,0].set_xlabel('Time (seconds)')
axs[0,0].set_ylabel('dF/F (%)')
axs[0,0].set_title('dF/F for rDLight')

lines = plot1#+ reward_ticks
labels = [l.get_label() for l in lines]  
legend = axs[0,0].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

axs[0,0].set_xlim(1000, 1060)
axs[0,0].set_ylim(-3, 7);

plot2=axs[1,0].plot(plotPoints, gACh_dF_F, 'g', label='gACh dF/F')
#

axs[1,0].set_xlabel('Time (seconds)')
axs[1,0].set_ylabel('dF/F (%)')
axs[1,0].set_title('dF/F for gACh')

lines = plot2#+ reward_ticks
labels = [l.get_label() for l in lines]  
legend = axs[1,0].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))

axs[1,0].set_xlim(1000, 1060)
ax[1,0].set_ylim(-10, 10);





################################################################            
    
            #subtract 405 from 465 and 506 traces (initial raw adjustment for movement/autofluorescence artifacts)
            isosFiltered_data465_LH = data465_LH-data405_LH
            isosFiltered_data560_LH = data560_LH-data405_LH
            isosFiltered_data465_RH = data465_RH-data405_RH
            isosFiltered_data560_RH = data560_RH-data405_RH
                    
            dataName_filteredFP = ['isosFiltered_data465_LH','isosFiltered_data560_LH','isosFiltered_data465_RH','isosFiltered_data560_RH']
            colors_filteredPlot = ['#008000','#FF0000','#90EE90','#FFC0CB']
            for f in range(len(dataName_filteredFP)):
                if eval(dataName_filteredFP[f]).any(): 
                    plt.figure(figsize=(10, 4))
                    plt.plot(eval(dataName_filteredFP[f])[500:], color=colors_filteredPlot[f])
                    plt.title('Time Series for {}'.format(dataName_filteredFP[f]))
                    plt.xlabel('approx. time (ms)')
                    plt.ylabel('raw fluorescence')
                    plt.show()    

#%% Part 1.3 - isolate epocs   
            epoc_ticker = dataTank.epocs.Tick.onset
            epoc_onsetCues = dataTank.epocs.sial.onset
            epoc_onsetITI = dataTank.epocs.sITI.onset
            epoc_beamBreak = dataTank.epocs.seam.onset

            # used Epoc marker as a duplicate (failed as of 5/1/24 - counting one second intervals)
            epoc_startTrial = dataTank.epocs.EpT_.onset
            epoc_startITI = dataTank.epocs.EpI_.onset
            epoc_startBeam = dataTank.epocs.EpC_.onset

# In[3]:

# Jupyter has a bug that requires import of matplotlib outside of cell with 
# matplotlib inline magic to properly apply rcParams
import matplotlib 
matplotlib.rcParams['font.size'] = 16 #set font size for all plots

# In[4]:
# ## Setup the variables for the data you want to extract
# 
# We will extract two different stream stores surrounding the 'PtAB' epoch event. We are interested in a specific event code for the shock onset.

REF_EPOC = 'PtAB' #event store name. This holds behavioral codes that are 
# read through ports A & B on the front of the RZ
SHOCK_CODE = [64959] #shock onset event code we are interested in

# make some variables up here to so if they change in new recordings you won't
# have to change everything downstream
ISOS = '_4054' # 405nm channel. Formally STREAM_STORE1 in maltab example
GCaMP = '_4654' # 465nm channel. Formally STREAM_STORE2 in maltab example
TRANGE = [-10, 20] # window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-10, -6] # baseline period within our window
ARTIFACT = float("inf") # optionally set an artifact rejection level

#call read block - new variable 'data' is the full data structure
data = tdt.read_block(BLOCKPATH)

# In[5]:
# ## Use epoc_filter to extract data around our epoc event
# 
# Using the 't' parameter extracts data only from the time range around our epoc event.<br>
# Use the 'values' parameter to specify allowed values of the REF_EPOC to extract.<br>
# For stream events, the chunks of data are stored in cell arrays structured as `data.streams[GCaMP].filtered`

data = tdt.epoc_filter(data, REF_EPOC, t=TRANGE, values=SHOCK_CODE)

# In[6]:

# Optionally remove artifacts. If any waveform is above ARTIFACT level, or
# below -ARTIFACT level, remove it from the data set.
total1 = np.size(data.streams[GCaMP].filtered)
total2 = np.size(data.streams[ISOS].filtered)

# List comprehension checking if any single array in 2D filtered array is > Artifact or < -Artifact
data.streams[GCaMP].filtered = [x for x in data.streams[GCaMP].filtered if not np.any(x > ARTIFACT) or np.any(x < -ARTIFACT)]

data.streams[ISOS].filtered = [x for x in data.streams[ISOS].filtered if not np.any(x > ARTIFACT) or np.any(x < -ARTIFACT)]

# Get the total number of rejected arrays
bad1 = total1 - np.size(data.streams[GCaMP].filtered)
bad2 = total2 - np.size(data.streams[ISOS].filtered)
total_artifacts = bad1 + bad2

# In[7]:
# Applying a time filter to a uniformly sampled signal means that the length of each segment could vary by one sample. Let's find the minimum length so we can trim the excess off before calculating the mean.

# More examples of list comprehensions
min1 = np.min([np.size(x) for x in data.streams[GCaMP].filtered])
min2 = np.min([np.size(x) for x in data.streams[ISOS].filtered])
data.streams[GCaMP].filtered = [x[1:min1] for x in data.streams[GCaMP].filtered]
data.streams[ISOS].filtered = [x[1:min2] for x in data.streams[ISOS].filtered]

# Downsample and average 10x via a moving window mean
N = 10 # Average every 10 samples into 1 value
F405 = []
F465 = []
for lst in data.streams[ISOS].filtered: 
    small_lst = []
    for i in range(0, min2, N):
        small_lst.append(np.mean(lst[i:i+N-1])) # This is the moving window mean
    F405.append(small_lst)

for lst in data.streams[GCaMP].filtered: 
    small_lst = []
    for i in range(0, min1, N):
        small_lst.append(np.mean(lst[i:i+N-1]))
    F465.append(small_lst)

#Create a mean signal, standard error of signal, and DC offset
meanF405 = np.mean(F405, axis=0)
stdF405 = np.std(F405, axis=0)/np.sqrt(len(data.streams[ISOS].filtered))
dcF405 = np.mean(meanF405)
meanF465 = np.mean(F465, axis=0)
stdF465 = np.std(F465, axis=0)/np.sqrt(len(data.streams[GCaMP].filtered))
dcF465 = np.mean(meanF465)

# In[8]:
# ## Plot epoc averaged response

# Create the time vector for each stream store
ts1 = TRANGE[0] + np.linspace(1, len(meanF465), len(meanF465))/data.streams[GCaMP].fs*N
ts2 = TRANGE[0] + np.linspace(1, len(meanF405), len(meanF405))/data.streams[ISOS].fs*N

# Subtract DC offset to get signals on top of one another
meanF405 = meanF405 - dcF405
meanF465 = meanF465 - dcF465

# Start making a figure with 4 subplots
# First plot is the 405 and 465 averaged signals
fig = plt.figure(figsize=(9, 14))
ax0 = fig.add_subplot(411) # work with axes and not current plot (plt.)

# Plotting the traces
p1, = ax0.plot(ts1, meanF465, linewidth=2, color='green', label='GCaMP')
p2, = ax0.plot(ts2, meanF405, linewidth=2, color='blueviolet', label='ISOS')

# Plotting standard error bands
p3 = ax0.fill_between(ts1, meanF465+stdF465, meanF465-stdF465,
                      facecolor='green', alpha=0.2)
p4 = ax0.fill_between(ts2, meanF405+stdF405, meanF405-stdF405,
                      facecolor='blueviolet', alpha=0.2)

# Plotting a line at t = 0
p5 = ax0.axvline(x=0, linewidth=3, color='slategray', label='Shock Onset')

# Finish up the plot
ax0.set_xlabel('Seconds')
ax0.set_ylabel('mV')
ax0.set_title('Foot Shock Response, %i Trials (%i Artifacts Removed)'
              % (len(data.streams[GCaMP].filtered), total_artifacts))
ax0.legend(handles=[p1, p2, p5], loc='upper right')
ax0.set_ylim(min(np.min(meanF465-stdF465), np.min(meanF405-stdF405)),
             max(np.max(meanF465+stdF465), np.max(meanF405+stdF405)))
ax0.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0]);

plt.close() # Jupyter cells will output any figure calls made, so if you don't want to see it just yet, close existing axis
            # https://stackoverflow.com/questions/18717877/prevent-plot-from-showing-in-jupyter-notebook
            # Note that this is not good code practice - Jupyter lends it self to these types of bad workarounds 


# ## Fitting 405 channel onto 465 channel to detrend signal bleaching
# 
# Scale and fit data. Algorithm sourced from Tom Davidson's Github:
# https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
# 

# In[9]:


Y_fit_all = []
Y_dF_all = []
for x, y in zip(F405, F465):
    x = np.array(x)
    y = np.array(y)
    bls = np.polyfit(x, y, 1)
    fit_line = np.multiply(bls[0], x) + bls[1]
    Y_fit_all.append(fit_line)
    Y_dF_all.append(y-fit_line)

# Getting the z-score and standard error
zall = []
for dF in Y_dF_all: 
   ind = np.where((np.array(ts2)<BASELINE_PER[1]) & (np.array(ts2)>BASELINE_PER[0]))
   zb = np.mean(dF[ind])
   zsd = np.std(dF[ind])
   zall.append((dF - zb)/zsd)
   
zerror = np.std(zall, axis=0)/np.sqrt(np.size(zall, axis=0))


# ## Heat Map based on z score of 405 fit subtracted 465

# In[10]:


ax1 = fig.add_subplot(412)
cs = ax1.imshow(zall, cmap=plt.cm.Greys, interpolation='none', aspect="auto",
                extent=[TRANGE[0], TRANGE[1]+TRANGE[0], 0, len(data.streams[GCaMP].filtered)])
cbar = fig.colorbar(cs, pad=0.01, fraction=0.02)

ax1.set_title('Individual z-Score Traces')
ax1.set_ylabel('Trials')
ax1.set_xlabel('Seconds from Shock Onset')

plt.close() # Suppress figure output again


# ## Plot the z-score trace for the 465 with std error bands

# In[11]:


ax2 = fig.add_subplot(413)
p6 = ax2.plot(ts2, np.mean(zall, axis=0), linewidth=2, color='green', label='GCaMP')
p7 = ax2.fill_between(ts1, np.mean(zall, axis=0)+zerror
                      ,np.mean(zall, axis=0)-zerror, facecolor='green', alpha=0.2)
p8 = ax2.axvline(x=0, linewidth=3, color='slategray', label='Shock Onset')
ax2.set_ylabel('z-Score')
ax2.set_xlabel('Seconds')
ax2.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0])
ax2.set_title('Foot Shock Response')

plt.close()


# ## Quantify changes as an area under the curve for cue (-5 sec) vs shock (0 sec)

# In[12]:


AUC = [] # cue, shock
ind1 = np.where((np.array(ts2)<-3) & (np.array(ts2)>-5))
AUC1= auc(ts2[ind1], np.mean(zall, axis=0)[ind1])
ind2 = np.where((np.array(ts2)>0) & (np.array(ts2)<2))
AUC2= auc(ts2[ind2], np.mean(zall, axis=0)[ind2])
AUC.append(AUC1)
AUC.append(AUC2)

# Run a two-sample T-test
t_stat,p_val = stats.ttest_ind(np.mean(zall, axis=0)[ind1],
                               np.mean(zall, axis=0)[ind2], equal_var=False)


# ## Make a bar plot

# In[13]:


ax3 = fig.add_subplot(414)
p9 = ax3.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)

# statistical annotation
x1, x2 = 0, 1 # columns indices for labels
y, h, col = max(AUC) + 2, 2, 'k'
ax3.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
p10 = ax3.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)

# Finish up the plot
ax3.set_ylim(0,y+2*h)
ax3.set_ylabel('AUC')
ax3.set_title('Cue vs Shock Response Changes')
ax3.set_xticks(np.arange(-1, len(AUC)+1))
ax3.set_xticklabels(['','Cue','Shock',''])

fig.tight_layout()
fig

