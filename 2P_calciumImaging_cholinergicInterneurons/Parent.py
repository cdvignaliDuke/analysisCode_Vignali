# -*- coding: utf-8 -*-
"""
Created on Tues Apr 23 11:35:25 2024

@author: Carlo Vignali
modified from 'Master.ipynp' written by Brandon Turner


"""

#%% Part 0.0 - Function import
# from suite2p.extraction import dcnv
import nbimporter
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
from PIL import ImageEnhance

# functions for processing data. 
from localFunctions_2PSLM import flatten
from localFunctions_2PSLM import _get_stim_list
from localFunctions_2PSLM import _cells_only
from localFunctions_2PSLM import _load_data
from localFunctions_2PSLM import _data_structure
from localFunctions_2PSLM import _norm_whole_trace
from localFunctions_2PSLM import _norm_by_stim

# from localFunctions_2PSLM import _suite2p_filter
from localFunctions_2PSLM import _rolling_filter
from localFunctions_2PSLM import _savgol_filter
from localFunctions_2PSLM import _reshape_df

from localFunctions_2PSLM import findspks
from localFunctions_2PSLM import find_peaks 
from localFunctions_2PSLM import _baseline_stats
from localFunctions_2PSLM import _pickle_waves
from localFunctions_2PSLM import round_to_nearest


#%% Part 0.1 - Metadata
#Define acquisition parameters
tau = .052 # timescale of indicator
framerate = 30 # CONSTANT fo frames per second. Based on period of 33.731ms. Do not change
frame_average = 1 # Rasters per frame. Adjust if you average frames *during acquisition*
fs = framerate/frame_average # sampling rate in Hz
neucoeff = 0.7 # neuropil coefficient (default is 0.7)
threshold = 3 # stdDev above baseline to find event

# Sample metadata
mouseID = 'ExCC106 1645910-NH' #Cage ID, CC#-(Toe ID)
mouse_geno = 'ChAT(+/cre) p'
group = 'pilot'
mouseDOB = '2023/12/02' # YYYY/MM/DD
mouseAGE = 110 #in DAYS
    
GCaMP_virus = 'IC, DIO-8f' #RO or IC. 
expression_time = 32 #in DAYS
    
wavelength = 920
objective = '20x'
sliceID = 'DS'
    
#set filename for output
experiment_date = '240315' #YYMMDD
experiment = 's4_picroACSF_BOT_10min'
    
#Parameters for analysis
stimulations = []
extra_pulses = 0 #write in list as ['x2', 'x3'], etc. Set to 0 if none. 
drug_treatments = ['PICRO'] #List drug treatments IN ORDER. ACSF is assumed for initial recordings. 
    
repetition   = 0 # number of times each stimulus is repeated
framenum     = 0 # number of frames for each stimulation epoch
stimtime     = 0 # time of stimulation after start of recording (ms)
    
#Constants for Smoothing Data
savgol_frames = 7
savgol_order = 3  
    
#Values for max/min suite2p filter. 
baseline = 'maximin' # take the running max of the running min after smoothing with gaussian
sig_baseline = 2 # in bins, standard deviation of gaussian with which to smooth (default 10)
win_baseline = 10 # in seconds, window in which to compute max/min filters (default 50)
    
ops = {'tau': tau, 'fs': fs, 'neucoeff': neucoeff, 'baseline': baseline, 'sig_baseline': sig_baseline, 'win_baseline': win_baseline}
doStim=int(input('Did stimulation or drug wash on occur in this recording?: '))
  
#In forSuite2P folder, list all experiments and then user inputs index for experiment of interest
currFolders = os.listdir("E:\\2PSLM\\CIN_gCaMP\\forSuite2P\\" + experiment_date); print(currFolders)
thisExp = currFolders[int(input('Which mouse from currFolders: '))-1]

#Set dataPath to suite2P output folder for that experiment
figPath = "E:\\2PSLM\\CIN_gCaMP\\forSuite2P\\" + experiment_date + '\\' + thisExp + '\\' + 'figures\\'
dataPath = "E:\\2PSLM\\CIN_gCaMP\\forSuite2P\\" + experiment_date + '\\' + thisExp + '\\' + 'suite2p\\' + 'plane0\\'
#%% Part 1.0 - Load data
#..............Previous iterations of the pipeline had been subtracting the flourescence of the neuropil (as calculated in suite2p) from the F signal before any processing. 
#..............As of 2021-08-24, removed this. Noticed it was yielding very low/negative values for some baselines and very high values for normalized dF/F (when F0 was <1 or <0, for example)
#..............To get the neuropil subtracted data, set 'subtract = True'

#CD to data directory
os.chdir(dataPath)
#Load data - see localFunctions_2PSLM for description of output variables
data, Fneu, iscell, stat, ops_orig, im, Fred, fchan = _load_data(subtract = False)
#
filelist = ops_orig['filelist']
print("Total Files: " + str(len(filelist)) + 'nOriginal Tiff Files in order')
for i in range(len(filelist)):
    print(filelist[i])


#%% Part 1.1 - Determine stimulation indices
#NOTE: This works if multi-pulse stimulation immediately follow single stim (i.e. 100, 100x2, 100x3). Else, if they repeat as such: [100, 100, 100, 100x2, 100x2, 100x2], 
# or it's a completely different order, you'll have to manually write it in manually for the 'stimulations' variable or change it later. 
##################################################################################################################################################################################################################
print(stimulations)
stimulations = _get_stim_list(stimulations, extra_pulses, drug_treatments)

##########################
# Manually written stims
##########################

# stimulations = ['900-pre','900-pre','900-pre','900-pre','900-pre',
#                'LED', 'LED-pre', 'LED-post', '900-inter',
#                'LED', 'LED-pre', 'LED-post', '900-inter',
#                'LED', 'LED-pre', 'LED-post', '900-inter',
#                'LED', 'LED-pre', 'LED-post', '900-inter',
#                'LED', 'LED-pre', 'LED-post', '900-inter',
#                '900-post', '900-post', '900-post', '900-post', '900-post', 
#                '900x3']

print(stimulations)

print(len(stimulations))

#%% Part 2.0 - remove non-cells from datalist
data_onlyCells, data_onlyCells_shape, indices = _cells_only(data, iscell)

#%% Part 2.1 - reformat datalist with stimulation protocol in consideration - irrelevant if no stimulation protocol or single stimulation protocol used
if doStim != 0:
    cellnum, stim_intensity, stim_num, baselineframenum = _data_structure(indices, stimulations, repetition, fs, stimtime=stimtime, chrimson=False)
else:
    cellnum = indices
    
#%% Part 3.0 - Data smoothing
row_size = len(indices)*len(stimulations)*repetition
print(row_size)
print(len(cellnum))
col_size = framenum

savgol_frames = 5 #Window for Savgol Polynomial Fit (7)
savgol_order = 2 #Order of Savgol Polynomial (2 or 3)
roll_window = 5 #Window for rolling mean calculations (7) 
roll_period = 1 #Min number of values required for average calculation. 

########################################
########## Smooth the data #############
########################################

wholeTrace_norm = _norm_whole_trace(d_cells) # whole trace normalized by trace mean
#e_norm = _norm_by_stim(d_cells, baselineframenum, row_size, col_size) # each epoch normalzized within itself
#e_norm_originalshape = np.reshape(e_norm, w_norm.shape) # whole trace normalized reshaped to match epoch format
#raw_e_reshape = np.reshape(d_cells, e_norm.shape) # Raw values reshaped to match epoch format

# sp = _suite2p_filter(e_norm, ops, adjust=False) # suite2p Min/Max filter of epoch format
wholeTrace_rollingFilter = _rolling_filter(w_norm, roll_window, roll_period) # Rolling mean window filter of each epoch
wholeTrace_savgolFilter = _savgol_filter(w_norm, savgol_frames, savgol_order) # Savgol filter of each epoch

#roll_raw = _rolling_filter(raw_e_reshape, roll_window, roll_period) # Rolling mean window filter of each epoch; Raw data, for QC. 
#sav_raw = _savgol_filter(raw_e_reshape, savgol_frames, savgol_order) # Savgol filter of each epoch; Raw data, for QC. 

sigma = 3
if doStim == 0:
    baselineframenum=0 # !! when no baseline period exists (i.e. spontaneous recording without manipulation)
else:
    baselineframenum = int(np.floor((stimtime-100)/1000*fs))

print(f"Baseline Frames: {baselineframenum}")

#
framesPerTick=50
roundedFrames = [round_to_nearest(len(data_onlyCells[0]))]
xVal = np.linspace(0,len(data_onlyCells[0]),len(data_onlyCells[0])+1)

for c in range(len(data_onlyCells)):
    fig,axes=plt.subplots(3,1,figsize=(12,21))
    fig.tight_layout(pad=3.0)

    axes[0].plot(xVal[1:],wholeTrace_norm[c],'k')
    axes[0].set_title('normalized trace - cell {}'.format(c))    
    axes[0].set_ylabel('raw fluorescence')
    
    axes[1].plot(xVal[1:],wholeTrace_rollingFilter[c],'k')
    axes[1].set_title('rolling filtered trace')    
    axes[1].set_ylabel('raw fluorescence')

    axes[2].plot(xVal[1:],wholeTrace_savgolFilter[c],'k')
    axes[2].set_title('savgol filtered trace')    
    axes[2].set_ylabel('raw fluorescence')

    plt.savefig('{directory}filterComparisons_{cell}.pdf'.format(directory=figPath, cell=str(c)), format='pdf', bbox_inches='tight', dpi=300)


#list(range(0,roundedFrames[0]+1,framesPerTick))
#%%

#.............In development, File '2021-05-20_s1_ChrimsonR_Sequence' yielding the following peaks. 
#.............Peaks from rolling == 1075
#.............Peaks from savgol  == 1161


    #Find spikes, Fmean, dFoverF using baseline normalized traces without filtering. 
    spkFrame_norm, spkTime_norm, spkBinary_norm, FmeanArr_norm, devFArr_norm = findspks(wholeTrace_norm, sigma, cellnum, baselineframenum, fs, 0)
    Fpeak_norm, dFoverF_norm, aucCalc_norm = find_peaks(wholeTrace_norm, baselineframenum, spkFrame_norm, spkTime_norm, spkBinary_norm, FmeanArr_norm, fs)
 
    #Find spikes, Fmean, dFoverF using baseline normalized traces passed through rolling filter. 
    spkFrame_rollingFilter, spkTime_rollingFilter, spkBinary_rollingFilter, FmeanArr_rollingFilter, devFArr_rollingFilter = findspks(wholeTrace_rollingFilter, sigma, cellnum, baselineframenum, fs, int(roll_window/2))
    Fpeak_rollingFilter, dFoverF_rollingFilter, aucCalc_rollingFilter = find_peaks(wholeTrace_rollingFilter, baselineframenum, spkFrame_rollingFilter, spkTime_rollingFilter, spkBinary_rollingFilter, FmeanArr_rollingFilter, fs)

    #Find spikes, Fmean, dFoverF using baseline normalized traces passed through savgol filter. 
    spkFrame_savgolFilter, spkTime_savgolFilter, spkBinary_savgolFilter, FmeanArr_savgolFilter, devFArr_savgolFilter = findspks(wholeTrace_savgolFilter, sigma, cellnum, baselineframenum, fs, savgol_frames)
    Fpeak_savgolFilter, dFoverF_savgolFilter, aucCalc_savgolFilter = find_peaks(wholeTrace_savgolFilter, baselineframenum, spkFrame_savgolFilter, spkTime_savgolFilter, spkBinary_savgolFilter, FmeanArr_savgolFilter, fs)
    
    #spk_frame2, spk_time2, spk_binary2, Fmean_arr2, devF_arr2 = findspks(roll_raw, sigma, cellnum, baselineframenum, fs, int(roll_window/2))
    #Fpeak2, dFoverF2, auc_calc2 = find_peaks(roll_raw, baselineframenum, spk_frame2, spk_time2, spk_binary2, Fmean_arr2, fs)
    #waves = roll_e
    
#..................................................................Q.C. -- Spike Count, Fpeak values, and dFoverF
print('Spikes from normalized traces: {}'.format(len(spk_binary[spk_binary == 1])))
print('Spikes from raw traces: {}'.format(len(spk_binary2[spk_binary2 == 1])))


#%%
waves_df = _pickle_waves(waves, iscell, stimulations, repetition, experiment, experiment_date)
waves_df.head(2)"


#%%
from localFunctions_2PSLM import _assign_id
from localFunctions_2PSLM import _cell_location

fredmean, frednorm, fmean_shaped, fnorm_shaped, d1, d1_strict, d1_d2, cutoff, interval_pos, interval_neg, raw_cutoff, raw_interval_pos, raw_interval_neg = _assign_id(Fred, iscell, stimulations, repetition)
xcoords_reshaped, ycoords_reshaped, im, xcoords, ycoords = _cell_location(stat, indices, stimulations, repetition, im)\n"


#%%
#Duplicate the metadata to fit into the dataframe
experiment_tag = np.repeat(str(experiment_date + '_'+ experiment),len(stimulations)*repetition*len(indices))
_date = np.repeat(experiment_date, len(stimulations)*repetition*len(indices))
_group = np.repeat(group, len(stimulations)*repetition*len(indices))

_mouseID = np.repeat(mouseID, len(stimulations)*repetition*len(indices))
_mouse_geno = np.repeat(mouse_geno, len(stimulations)*repetition*len(indices))
_mouseDOB = np.repeat(mouseDOB, len(stimulations)*repetition*len(indices))
_mouseAGE = np.repeat(mouseAGE, len(stimulations)*repetition*len(indices))
_GCaMP_virus = np.repeat(GCaMP_virus, len(stimulations)*repetition*len(indices))
_expression_time = np.repeat(expression_time, len(stimulations)*repetition*len(indices))
_wavelength = np.repeat(wavelength, len(stimulations)*repetition*len(indices))
_objective = np.repeat(objective, len(stimulations)*repetition*len(indices))
_sliceID = np.repeat(sliceID, len(stimulations)*repetition*len(indices))


d = {'cell_num': cellnum, 
     'Experiment': experiment_tag,
     'Date': _date,
     'Group': _group,
     'ID':d1, 
     'StrictID': d1_strict,
     'Red Intensity (Avg)': fmean_shaped, 
     'Red norm': fnorm_shaped,
     'x-coords': xcoords_reshaped, 
     'y-coords': ycoords_reshaped, 
     'stim_intensity': stim_intensity, 
     'stim_num': stim_num, 
     'F_mean': Fmean_arr, 
     'StdDev' : devF_arr, 
     'spike_happened': spk_binary, 
     'spike_time(ms)': spk_time, 
     'dF': Fpeak, 
     'AUC': auc_calc, 
     'dF/F(%)': dFoverF,
     'MouseID': _mouseID,
     'Genotype': _mouse_geno, 
     'DOB': _mouseDOB, 
     'Age (days)': _mouseAGE, 
     'GCaMP Virus': _GCaMP_virus, 
     'Expression Time (days)': _expression_time, 
     'Excitation Wavelength (nm)': _wavelength, 
     'Objective': _objective, 
     'Slice ID': _sliceID
    }

df = pd.DataFrame(data=d)
df.insert(6, 'Stim Time', int(stimtime))
df['ID'] = np.where((df.ID==1), 'D1', 'D2')


export_name = experiment_date + '_' + experiment + '.csv'
export_csv = df.to_csv (export_name, index = None, header=True)
print(export_name)

df[df['spike_happened'] == 1].head(3)"


#%%
from localFunctions_2PSLM import _plot_fov
from localFunctions_2PSLM import _plot_id

_plot_id(frednorm, fredmean, d1_d2, cutoff, interval_pos, interval_neg, raw_cutoff, raw_interval_pos, raw_interval_neg, experiment, experiment_date)
_plot_fov(fchan, ops_orig, im, indices, experiment_date, experiment, xcoords, ycoords, d1_d2)"


#%%
from localFunctions_2PSLM import _plot_normed_traces
from localFunctions_2PSLM import _plot_filtered_traces
from localFunctions_2PSLM import _plot_savgol_windows
from localFunctions_2PSLM import _plot_savgol_order
from localFunctions_2PSLM import _norm_raw_corr
from localFunctions_2PSLM import _plot_found_spikes       
from localFunctions_2PSLM import _wave_plot

#..........................................................................................................................................................

roi_num = 15
_plot_normed_traces(w_norm, e_norm_originalshape, d_cells, roi_num, framenum, stimulations, repetition, savefig=True) 
#...................................Plot normalized traces (whole trace norm and by_epoch norm)
# _plot_savgol_windows(e_norm, roi_num, framenum, d_cells.shape, stimulations, repetition, xmin = 1200, xlim = 1500) 
#....................................default order is 3. Windows are 7, 15, & 29. 
# _plot_savgol_order(e_norm, roi_num, framenum, d_cells.shape, stimulations, repetition, xmin = 1200, xlim = 1500) 
#......................................default window is 7. Order is 2 & 3. 
sp = roll_e.copy()
_plot_filtered_traces(sp, roll_e, sav_e, roi_num, framenum, d_cells.shape, stimulations, repetition, xmin=None, xlim=None) 
#..........................................compare suite2p min/max, rolling average, and Savgol filters. 


#..........................................................................................................................................................

cells_plotted = [i for i in range(15,16)]
_plot_found_spikes(roll_e, sav_e, d_cells, df, cells_plotted, spk_frame, Fmean_arr, Fpeak, 
#...............................................................Plots traces from specified ROIs with identified spikes. 
                   indices, stimulations, repetition, fs, stimtime, framenum, framerate,
                   rolling, 
                   savefig = False) 

#..........................................................................................................................................................
if rolling:
    title = \"Rolling\"
else: 
    title = \"Savgol\"
_norm_raw_corr(Fpeak, Fpeak2, dFoverF, dFoverF2, savefig=True, title=title) 
#.............................................................................Plot comparison of peak amplitudes from Raw and Normalized Values

#..........................................................................................................................................................

#................................For Chrimson Experiments. 
# prestim = [1,2,3,4,5]
# led = [6 + 4*i for i in range(5)]
# pre = [7 + 4*i for i in range(5)]
# post = [8 + 4*i for i in range(5)]
# stim = [9 + 4*i for i in range(5)]
# poststim = [26, 27, 28, 29, 30]
# check = [31]
#................................

first_roi = 15
last_roi = 16

stim = [1, 2, 3, 4, 5, 11, 12, 13, 14, 15]

_wave_plot(waves_df, stim, first_roi, last_roi) 
#..........................................................................................................Plots Waveform from specified ROIs ('first_roi', 'last_roi')"


#%%


#%%


#%%


#%%


#%%


#%%


#%%


