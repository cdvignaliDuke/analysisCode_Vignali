# -*- coding: utf-8 -*-
"""
Created on Tues Apr 23 11:35:25 2024

@author: Carlo Vignali
modified from 'Functions.ipynp' written by Brandon Turner


"""

#%% 
#from suite2p.extraction import dcnv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
from PIL import ImageEnhance

#%%
def flatten(l):
    flatList = []
    for elem in l:
    # if an element of a list is a list
    # iterate over this list and add elements to flatList
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

#%%
def _get_stim_list(pulses, extra_pulses, drug_treatments):
    stimulations = pulses.copy()
    if extra_pulses != 0: 
        pulses = []
        for x in range(len(extra_pulses)):
            pulses.append(x)    
        for y in extra_pulses:
            pulses[extra_pulses.index(y)] = [str(i) + y for i in stimulations]
    else:
        pulses = []

    stims_acsf = []
#stimulations = stims_acsf
    for i in range(len(stimulations)):
        stims_acsf.append(stimulations[i])
        for j in range(len(pulses)):
            stims_acsf.append(pulses[j][i])

    stims_acsf = flatten(stims_acsf)

    stim_drugs = []
    for x in drug_treatments:                
        stim_drug_temp = [(str(i)+'-' + x) for i in stims_acsf]
        stim_drugs.append(stim_drug_temp)

    stimulations = stims_acsf
    stimulations.extend(stim_drugs)
    stimulations = flatten(stimulations)
    print("Total Stims: {len(stimulations)} + /nStim Vals:/n {stimulations}")
    return stimulations

#%%
def _load_data(coeff = 0.7, subtract = False): 
    neucoeff = coeff
    print("Returns values: 'data', 'Fneu', 'iscell', 'stat', 'Fc', 'ops_orig', 'im'")
    F = np.load('F.npy')           #raw fluorescence traces
    Fneu = np.load('Fneu.npy')     #neuropil fluorescence
    iscell = np.load('iscell.npy') #cell classifier
    stat = np.load('stat.npy', allow_pickle=True)     #aggregate statistics for each ROI
    ops_orig = np.load('ops.npy', allow_pickle=True).item()
    im = np.zeros((ops_orig['Ly'], ops_orig['Lx']))
    if subtract: 
        Fc = F - neucoeff * Fneu       #subtract neuropil component from each ROI 
        data = Fc
    else: 
        data = F
        
    try: 
        Fred = np.load('F_chan2.npy')
        fchan = True
        print("Data Shape: {data.shape}")
        return data, Fneu, iscell, stat, ops_orig, im, Fred, fchan

    except: 
        fchan = False
        Fred = np.empty(F.shape)
        print("No red channel")
        print("Data Shape: {data.shape}")
        return data, Fneu, iscell, stat, ops_orig, im, Fred, fchan


#%%
def _cells_only(data, iscell):
    indices, prob = np.where(iscell==1)
#     print(indices)
    d = data.copy()
    print("Original Shape: {d.shape}")
    d = d[iscell[:,0]==1]
    print("New Shape: {d.shape}")
    trimmed_shape = d.shape
    return d, trimmed_shape, indices

#%%
def _norm_whole_trace(data): # Normalize data as F/F0. 
    norm_trace = np.zeros(data.shape)
    for i in range(len(data)):
        for x in range(len(data[i,:])):
            norm_trace[i,x] = (data[i,x])/np.mean(data[i,:])
    return norm_trace


#%%
def _norm_by_stim(data, baselineframenum, rows, columns):
    print("Returns reshaped dataframe normalized within each stim period. Must specify 'data', 'baselineframenum', 'rows', and 'columns'")
    d = data.copy()
    d = np.reshape(d, (rows, columns))
    norm_trace = np.zeros(d.shape)
    for i in range(len(d)):
        for x in range(len(d[i,:])):
            norm_trace[i,x] = (d[i,x])/np.mean(d[i,:baselineframenum])
    return norm_trace  

#%%
def _suite2p_filter(a, ops, adjust=True):
    print("Takes input 'd' (data) & 'ops'")
    d = a.copy()
    normt_filtered = dcnv.preprocess(
         F=d,
         baseline=ops['baseline'],
         win_baseline=ops['win_baseline'],
         sig_baseline=ops['sig_baseline'],
         fs=ops['fs'],
         prctile_baseline=100
     )
    if adjust: #..................................adjust so that baseline is at 1, not zero. 
        adjustment = np.ones(normt_filtered.shape)            
        normt_filtered = np.add(normt_filtered, adjustment)
    print("Output Shape: {normt_filtered.shape}")
    return normt_filtered

#%%
def _rolling_filter(a, window, min_vals):
    b = a.copy()
    c = pd.DataFrame(data=b).astype(float)
    d = c.rolling(window, min_periods = min_vals, axis = 1).mean()
    e = d.to_numpy()
    print("Returned Value Shape: {e.shape}")
    return e

#%% findpks - Original spike finder. Gets 1st point to cross threshold for each stim epoch
def findspks (traces, sigma, cellnum, baselineframenum, fs, savgol_framenum):# 
    spike_frame = np.full(len(cellnum) ,np.nan)  # frame number of the first spike detected
    Fmean_arr = np.full(len(cellnum) ,np.nan)  # frame number of the first spike detected
    devF_arr = np.full(len(cellnum) ,np.nan)  # frame number of the first spike detected
    sigma = sigma
    for i in range(len(cellnum)):
        data = traces[i, baselineframenum::]
        Fmean = np.mean(traces[i,savgol_framenum:baselineframenum-savgol_framenum])
        devF = np.std(traces[i,savgol_framenum:baselineframenum-savgol_framenum])
        
        Fmean_arr[i] = Fmean
        devF_arr[i] = devF
        
        thr = 1.2*Fmean #finding spikes with 20 percent more than Fmean
        std_thr = Fmean + sigma * devF                     
        data_thr = np.where((data>thr) & (data>std_thr))[0]
        if data_thr.size>0: 
            spike_frame[i] = int(data_thr[0]+baselineframenum)
            
    spike_time = spike_frame*1000/fs  #-stimtime    # converting to miliseconds from the stimulus
    spike_time = np.where(spike_time<0, np.nan, spike_time) # removing negative timings
    spike_binary = ~np.isnan(spike_time)*1
    return spike_frame, spike_time, spike_binary, Fmean_arr, devF_arr

#%%
def find_peaks (data, baselineframenum, spike_frame, spike_time, spike_binary, Fmean, fs, search_window = 500, plot_floor = False, trimmed_shape = (0,0)):
    d = data.copy()
    Fpeak = np.empty(spike_frame.shape)
    for i in range(len(spike_frame)):
        if ~np.isnan(spike_frame[i]):
            Fpeak[i] = np.max(d[i,int(spike_frame[i]):int(spike_frame[i] + (search_window/1000*fs))])

    Fpeak = np.where(spike_binary==1, Fpeak, np.nan)
    dFoverF = (Fpeak-Fmean)/Fmean * 100

    adjustments = np.ones(d.shape)
    floor = np.subtract(d, adjustments)
    auc_calc = np.sum(floor[:, baselineframenum:int(np.floor(search_window*2/1000*fs))], axis=1)  # Calculated AUC in each sweep. Up to 2x peak search window. 
    auc_calc= np.where(spike_binary==1, auc_calc, np.nan)
    
    print("Spike array shape: {spike_binary.shape}")
    print("FPeak array shape: {Fpeak.shape}")
    
    if plot_floor: #................................................Floored traces to be plotted. QC  use only.  
    #4) Find area under the curve (AUC)
        adjustments = np.ones(baselined.shape)
        floor = np.subtract(baselined, adjustments)
        auc = np.where(floor<0, 0, floor)     
        
        a_floor = np.reshape(auc, trimmed_shape)
        
        fig = plt.figure(dpi=150)
        n = rand.randint(range(len(a_floor)))
        plt.plot(a_floor[n], linestyle = 'dotted', marker = 'o', color = 'tab:gray', alpha = 0.7)
    return Fpeak, dFoverF, auc_calc

#%%
def _data_structure(indices, stimulations, repetition, fs, chrimson = True, stimtime=1000):#.................Build data structure
    cellnum =  np.repeat(indices, len(stimulations)*repetition)
    stim_intensity = np.tile(np.repeat(stimulations, repetition), len(indices))
    if chrimson:
        chr_stims = [1,2,3,4,5,
                     1,1,1,1,
                     2,2,2,2,
                     3,3,3,3,
                     4,4,4,4,
                     5,5,5,5,
                     1,2,3,4,5,
                     1]
        print('Chrimson Experiment')
        print('Stims: {chr_stims}')
        stim_num = np.tile(chr_stims, len(indices))
    else:
        stim_num = np.tile(range(repetition), len(indices)*len(stimulations))
        stim_num = stim_num +1

    baselineframenum = int(np.floor((stimtime-100)/1000*fs))
    print(stimulations)
    return cellnum, stim_intensity, stim_num, baselineframenum

#%%
def round_to_nearest(x, bases=(50, 100)):
    return min((base * round(float(x)/base) for base in bases), key=lambda v: abs(v-x))

#%%
def _savgol_filter(data, savgol_frames, savgol_order):#.......................Run Savgol (polynomial) filtering on data
    f = savgol_filter(data, savgol_frames,savgol_order)
    return f

#%%
def _reshape_df(data, cellnum, framenum):
    reshaped = np.reshape(data, (len(cellnum), framenum))
    return reshaped

#%%
#3) Find Fmean in each sweep/stim
def _baseline_stats(data, baselineframenum):
    Fmean = np.mean(baselinedreshaped[:,0:baselineframenum], axis=1)
    devF = np.std(baselinedreshaped[:, 0:baselineframenum], axis=1)
    return Fmean, defF

#%%
def _pickle_waves(data, iscell, stimulations, repetition, experiment, experiment_date):
    d = data.copy()

    ic = iscell[:,0].copy()
    for i in range(len(ic)):
        if ic[i] >= 1: 
            ic[i] = int(i)
        else: 
            ic[i] = np.nan
    ic = ic[~np.isnan(ic)]
    
    colnums = len(stimulations) * repetition
    rownums = len(ic)
    xvals = range(len(d[1]))
    
    val_array = np.empty(shape = (rownums, colnums+1), dtype = object)
    for n in range(rownums): 
        t1 = d[n*colnums:(n+1)*colnums,:]
        for s in range(len(t1[:,0])):
            t2 = t1[s,:]
            val_array[n, s+1] = [xvals, t2]
    val_array[:,0] = np.nan
    
    colheader = ["ROI#"]
    for i in range(len(stimulations)*repetition):
        colheader.append("Stim"+str(i+1))

    df_traces = pd.DataFrame(data = val_array, columns = colheader)
    df_traces['ROI#'] = ["{experiment_date}_{experiment}-Cell_Num-{int(ic[i])}" for i in range(len(ic))]
#     '2021-05-20_s1_ChrimsonR_Sequence-Cell_Num-10'
    df_traces.to_pickle("{experiment_date}_{experiment}_waveforms.zip", compression = "infer")
    return df_traces

#%%
def _plot_filtered_traces(suite2p_arr, roll_arr, sav_arr, k, framenum, shape, stimulations, repetition, xmin = None, xlim = None, savefig=True):
    a = np.reshape(suite2p_arr.copy(), shape)
    b = np.reshape(roll_arr.copy(), shape)
    c = np.reshape(sav_arr.copy(), shape)
   
    fig = plt.figure(dpi=150, figsize = (10,5))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax1.plot(a[k], alpha = 0.5, color = 'k')
    ax2.plot(b[k], alpha = 0.5, color = 'b')
    ax3.plot(c[k], color = 'k', alpha = 0.5)
    for i in range(len(stimulations)*repetition):
        x = framenum + framenum*i
        if i%2 == 0:
            ax1.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax3.axvline(x, color='r', zorder = 0, alpha = 0.2)
        else: 
            ax1.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax3.axvline(x, color='k', zorder = 0, alpha = 0.2)
        
    axmax = np.max([ax1.get_ylim(), ax2.get_ylim(), ax3.get_ylim()])
    axmin = np.min([ax1.get_ylim(), ax2.get_ylim(), ax3.get_ylim()])
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)
    ax3.set_ylim(axmin, axmax)
    
    ax1.set_ylabel('F/F0')
    ax1.axhline(y=0, color = 'b', zorder = 0, alpha = 0.5)
    ax2.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)
    ax3.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)

    ax1.set_title("MaxMin Filter")
    ax2.set_title("Rolling Mean")
    ax3.set_title("Savgol Filter")
    ax3.set_ylabel("F (A.U.)")
    
    ax1.set_xlim(xmin, xlim)
    ax2.set_xlim(xmin, xlim)
    ax3.set_xlim(xmin, xlim)
    
    plt.tight_layout()
    if savefig:
        plt.savefig("FilterComparison.png", bbox_inches='tight')

#%%
def _plot_savgol_windows(raw, k, framenum, shape, stimulations, repetition, window1=7, window2=15, window3=29, order=3, xmin = None, xlim = None, savefig=True):
    a = _savgol_filter(raw, window1, order)
    b = _savgol_filter(raw, window2, order)
    c= _savgol_filter(raw, window3, order)
    
    a = np.reshape(a, shape)
    b = np.reshape(b, shape)
    c = np.reshape(c, shape)
    d = np.reshape(raw.copy(), shape)
    
    fig = plt.figure(dpi=150, figsize = (10,5))
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)
    
    ax1.plot(a[k], alpha = 0.5, color = 'k')
    ax2.plot(b[k], alpha = 0.5, color = 'b')
    ax3.plot(c[k], color = 'k', alpha = 0.5)
    ax4.plot(d[k], color = 'k', alpha = 1)
    for i in range(len(stimulations)*repetition):
        x = framenum + framenum*i
        if i%2 == 0:
            ax1.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax3.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax4.axvline(x, color='r', zorder = 0, alpha = 0.2)
        else: 
            ax1.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax3.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax4.axvline(x, color='k', zorder = 0, alpha = 0.2)
        
    axmax = np.max([ax1.get_ylim(), ax2.get_ylim(), ax3.get_ylim(), ax4.get_ylim()])
    axmin = np.min([ax1.get_ylim(), ax2.get_ylim(), ax3.get_ylim(), ax4.get_ylim()])
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)
    ax3.set_ylim(axmin, axmax)
    ax4.set_ylim(axmin, axmax)
    
    ax1.set_xlim(xmin, xlim)
    ax2.set_xlim(xmin, xlim)
    ax3.set_xlim(xmin, xlim)
    ax4.set_xlim(xmin, xlim)

    ax1.axhline(y=0, color = 'b', zorder = 0, alpha = 0.5)
    ax2.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)
    ax3.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)
    ax4.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)

    ax1.set_title("Savgol Filter: 7,3")
    ax2.set_title("Savgol Filter: 15,3")
    ax3.set_title("Savgol Filter: 30,3")
    ax4.set_title("Normalized; Unfiltered")
    
    ax1.set_ylabel('F/F0')
    ax2.set_ylabel('F/F0')
    ax3.set_ylabel('F/F0')
    ax4.set_ylabel('F/F0')
    
    plt.tight_layout()
    if savefig:
        plt.savefig("Savgol_Filter_WindowComparison.png", bbox_inches='tight')


#%%
def _plot_savgol_order(raw, k, framenum, shape, stimulations, repetition, window = 7, order1 = 2, order2 = 3,  xmin = None, xlim = None, savefig=True):
    
    a = _savgol_filter(raw, window, order1)
    b = _savgol_filter(raw, window, order2)
    
    a = np.reshape(a, shape)
    b = np.reshape(b, shape)
    d = np.reshape(raw.copy(), shape)  
    
    fig = plt.figure(dpi=150, figsize = (10,5))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
#     ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(313)
    
    ax1.plot(a[k], alpha = 0.5, color = 'k')
    ax2.plot(b[k], alpha = 0.5, color = 'b')
#     ax3.plot(c[k], color = 'k', alpha = 0.5)
    ax4.plot(d[k], color = 'k', alpha = 1)
    for i in range(len(stimulations)*repetition):
        x = framenum + framenum*i
        if i%2 == 0:
            ax1.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='r', zorder = 0, alpha = 0.2)
#             ax3.axvline(x, color='r', zorder = 0, alpha = 0.2)
            ax4.axvline(x, color='r', zorder = 0, alpha = 0.2)
        else: 
            ax1.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax2.axvline(x, color='k', zorder = 0, alpha = 0.2)
#             ax3.axvline(x, color='k', zorder = 0, alpha = 0.2)
            ax4.axvline(x, color='k', zorder = 0, alpha = 0.2)
        
    axmax = np.max([ax1.get_ylim(), ax2.get_ylim(),  ax4.get_ylim()])
    axmin = np.min([ax1.get_ylim(), ax2.get_ylim(),  ax4.get_ylim()])
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)
#     ax3.set_ylim(axmin, axmax)
    ax4.set_ylim(axmin, axmax)
    
    ax1.set_xlim(xmin, xlim)
    ax2.set_xlim(xmin, xlim)
#     ax3.set_xlim(xmin, xlim)
    ax4.set_xlim(xmin, xlim)

    ax1.axhline(y=0, color = 'b', zorder = 0, alpha = 0.5)
    ax2.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)
#     ax3.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)
    ax4.axhline(y=1, color = 'b', zorder = 0, alpha = 0.5)

    ax1.set_title("Savgol Filter: 7,2")
    ax2.set_title("Savgol Filter: 7,3")
#     ax3.set_title("Savgol Filter: 30,2")
    ax4.set_title("Normalized; Unfiltered")
    
    ax1.set_ylabel('F/F0')
    ax2.set_ylabel('F/F0')
#     ax3.set_ylabel('F/F0')
    ax4.set_ylabel('F/F0')
    
    plt.tight_layout()
    if savefig:
        plt.savefig("Savgol_Filter_OrderComparison.png", bbox_inches='tight')


#%%
def _plot_normed_traces(by_cell, by_epoch, raw, k, framenum, stimulations, repetition, savefig=True):
    a = by_cell.copy()
    b = by_epoch.copy()
    c = raw.copy()
    
    fig = plt.figure(dpi=150)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(212)
    ax1.plot(a[k], alpha = 0.5, color = 'k')
    ax2.plot(b[k], alpha = 0.5, color = 'b')
    ax3.plot(c[k], color = 'k', alpha = 0.5)
    for i in range(len(stimulations)*repetition):
        x = framenum + framenum*i
        ax1.axvline(x, color='r', zorder = 0, alpha = 0.2)
        ax2.axvline(x, color='r', zorder = 0, alpha = 0.2)
        ax3.axvline(x, color='r', zorder = 0, alpha = 0.2)
    axmax = np.max([ax1.get_ylim(), ax2.get_ylim()])
    axmin = np.min([ax1.get_ylim(), ax2.get_ylim()])
    ax1.set_ylim(axmin, axmax)
    ax2.set_ylim(axmin, axmax)
    ax1.set_ylabel('F/F0')
    ax1.axhline(y=0, color = 'b', zorder = 0, alpha = 0.5)
    ax2.axhline(y=0, color = 'b', zorder = 0, alpha = 0.5)

    ax1.set_title("Norm whole trace")
    ax2.set_title("Norm each epoch")
    ax3.set_title("Raw Flourescence values")
    ax3.set_ylabel("F (A.U.)")
    plt.tight_layout()
    if savefig:
        plt.savefig("NormedTraceComparison.png", bbox_inches='tight')


#%%
def _norm_raw_corr (Fpeak1, Fpeak2, dFoverF1, dFoverF2, savefig = True, title="Rolling"): 
    import seaborn as sns
    Fpeak = Fpeak1.copy()
    Fpeak2 = Fpeak2.copy()
    dFoverF = dFoverF1.copy()
    dFoverF2 = dFoverF2.copy()
    
    fig = plt.figure(dpi = 150, figsize = (8,4))
    grid = plt.GridSpec(2, 10, wspace=5, hspace=.5)

    ax1 = plt.subplot(grid[0, 0:2])
    ax2 = plt.subplot(grid[1, 0:2])
    ax3 = plt.subplot(grid[:, 3:6])
    ax4 = plt.subplot(grid[:, 6:])

    ax1.plot(np.sort(Fpeak))
    ax1.set_title('Norm Peaks')
    ax2.plot(np.sort(Fpeak2))
    ax2.set_title('Raw Peaks')
    ax3.scatter(dFoverF, dFoverF2)
    x = np.linspace(0, np.max(dFoverF[~np.isnan(dFoverF)]), 20)

    ax3.plot(x,x, linestyle = 'dotted', c='r')
    ax3.set_xlabel('dF/F(%) Normed')
    ax3.set_ylabel('dF/F(%) Raw')
#     ax3.set_xticks([0, 500, 1000, 1500, 2000])
#     ax3.set_xticklabels([0, 500, 1000, 1500, 2000], rotation = 30)
#     ax3.set_yticks([0, 500, 1000, 1500, 2000])
#     ax3.set_yticklabels([0, 500, 1000, 1500, 2000], rotation = 30)
    ax3.set_title("1:1 Correlation Expected")
    lim_max = np.max([ax3.get_ylim(), ax3.get_xlim()])

    sns.histplot(dFoverF2, bins = 150, ax = ax4)
    ax4.set_xlim(0, 600)
    ax4.set_xlabel("dF/F(%) Raw Peaks")
    ax4.set_title("Histogram; Bins = 300")
    if savefig: 
        plt.savefig("SpikePeak_QC_"+title+".png", bbox_inches='tight')


#%%
def _wave_plot(waves_df, stimplot, roi_start, roi_end, interval=1):
    df_plot = waves_df.copy()
    rois = df_plot['ROI#'].drop_duplicates().to_list()
    for i in rois[roi_start:roi_end]:
        plotdf = df_plot[df_plot['ROI#'] == i]
        fig = plt.figure()
        for n in stimplot:
            pdf = plotdf.iloc[:,n]
            plt.plot(pdf[rois.index(i)][0], pdf[rois.index(i)][1], label = n, c = 'g', alpha = 0.5)
            plt.suptitle("ROI#"+str(i))
            plt.legend()


#%%
def _plot_fov(fchan, ops_orig, im, indices, experiment_date, experiment, xcoords, ycoords, d1_d2):
    x = xcoords.copy()
    y = ycoords.copy()
    ncells = indices
    #Load image of red channel
    if fchan: 
        mean_im = ops_orig['meanImg_chan2']
    else: 
        mean_im = ops_orig['meanImg']

    #Create a plot
    fig1 = plt.figure(figsize  = (7,7), dpi = 150)
    plt.suptitle(experiment_date + '--' + experiment)

    rois_only = fig1.add_subplot(221)
    redimage = fig1.add_subplot(222)
    ax_image = fig1.add_subplot(223)
    ax = fig1.add_subplot(224)

    #Plot the ROIs and overlay images
    rois_only.imshow(im, cmap = 'Blues', vmax = 1)
    rois_only.set_title('ROIs')

    #Show mean red image (change vmax if too washed out)
    if fchan: 
        redimage.imshow(mean_im, cmap='Reds', vmax=1000)
    else: 
        redimage.imshow(mean_im, cmap='Greens', vmax=np.max(mean_im)/5)

    # Show overlay w/ masked image
    mask = np.zeros(shape = im.shape)
    mask = np.where(im > 0, 1, 0)
    masked = im*mask

    # ax_rois.imshow(masked, cmap='magma', vmax=1, vmin=0, alpha = 0.7)
    ax_image.set_xlabel('Location (pixels)', color='black')
    ax_image.set_ylabel('Location (pixels)', color='black')
    if fchan: 
        ax_image.imshow(mask*mean_im, cmap='Reds', vmin=0, vmax=1000)
        ax_image.set_title('Red overlay')
    else: 
        ax_image.imshow(mask*mean_im, cmap='Greens', vmin=0, vmax=np.max(mean_im)/5)
        ax_image.set_title('GCaMP ONLY; no red channel')

    ax.scatter(x, y, s=20, c=d1_d2, cmap='prism_r' )
    ax.set_xlabel('Location (pixels)', color='black')
    ax.set_ylabel('Location (pixels)', color='black')
    ax.set_ylim(ymin=0, ymax=512)
    ax.set_xlim(0, 512)
    ax.invert_yaxis()

    #Annotate graphs with ROI/cell numbers 
    annotation = zip(ncells, x, y)
    for n, x, y in annotation:
        ax.annotate(n, xy=(x,y), color = 'black')

    plt.savefig(experiment_date + '--' + experiment + '--ROIs.png', bbox_inches='tight')
    plt.show()


#%%
def _cell_location(stat, indices, stimulations, repetition, im):
    ncells = indices.copy()
    centers = []
    for n in ncells:
        centers += [(stat[n]['med'])]
    centers = np.asarray(centers)
    xcoords = centers[:,1]
    ycoords = centers[:,0]
    x = xcoords.tolist()
    y = ycoords.tolist()
    print(x)
    print(y)
    #Identify centroids of selected cells (above)
    xcoords_reshaped = np.repeat(xcoords, len(stimulations)*repetition)
    ycoords_reshaped = np.repeat(ycoords, len(stimulations)*repetition)
    xcoords_reshaped.shape

    ypix_mid = []
    xpix_mid = []
    for n in ncells:
        ypix = stat[n]['ypix'][~stat[n]['overlap']]
        xpix = stat[n]['xpix'][~stat[n]['overlap']]
        #ypix_mid.append(sum(ypix)/(1 if len(ypix)==0 else len(ypix)))
        #xpix_mid.append(sum(xpix)/(1 if len(xpix)==0 else len(xpix)))  
        im[ypix,xpix] = n+1
    return xcoords_reshaped, ycoords_reshaped, im, x, y


#%%
def _assign_id(Fred, iscell, stimulations, repetition):
#Extract the data from the red channel (channel 2)
    indices, prob = np.where(iscell==1)
    ncells = indices.copy()
    fredmean = np.mean(Fred, axis = 1)
    
    print("F-RED array shape: ", Fred.shape)
    print("F-RED mean array shape: ", fredmean.shape)
    fredmean = fredmean[iscell[:,0]==1]
    frednorm = [i/np.mean(fredmean, axis = 0) for i in fredmean]

    #Reshape Arrays
    fmean_shaped = np.repeat(fredmean, len(stimulations)*repetition)
    fnorm_shaped = np.repeat(frednorm, len(stimulations)*repetition)


    #Define a cutoff (as a range on the minimum and maximum red fluorescence value) to define D1 vs D2
    #cutoff is the median. Split cells by confidence interval. 
    cutoff = np.median(frednorm, axis = 0)
    len_cells = len(ncells)

    interval_pos = (1 + len_cells/2) + (1.96*np.sqrt(len_cells)/2)
    interval_pos = np.sort(frednorm)[int(interval_pos.round())]

    interval_neg = (len_cells/2)-(1.96*np.sqrt(len_cells)/2)
    interval_neg = np.sort(frednorm)[int(interval_neg.round())]

    #similar cutoffs for raw values
    raw_cutoff = np.median(fredmean, axis=0)
    raw_interval_pos = (1 + len_cells/2) + (1.96*np.sqrt(len_cells)/2)
    raw_interval_pos = np.sort(fredmean)[int(raw_interval_pos.round())]

    raw_interval_neg = (len_cells/2)-(1.96*np.sqrt(len_cells)/2)
    raw_interval_neg = np.sort(fredmean)[int(raw_interval_neg.round())]


    #store D1/D2 ID in an array
    d1 = [1 if i> cutoff else 0 for i in fnorm_shaped]  #.................................................................................................for strict ID purposes
    d1_strict = ['D1' if i>=(interval_pos) else 'D1?' if interval_pos>i>cutoff else 'D2?' if cutoff>i>interval_neg else 'D2' for i in fnorm_shaped] #.....for strict ID purposes
    d1_d2 = [2 if i>=interval_pos else 1 if interval_pos>i>interval_neg else 0 for i in frednorm]  #......................................................for coloring in graphs
    
    return fredmean, frednorm, fmean_shaped, fnorm_shaped, d1, d1_strict, d1_d2, cutoff, interval_pos, interval_neg, raw_cutoff, raw_interval_pos, raw_interval_neg


#%%
def _plot_id(frednorm, fredmean, d1_d2, cutoff, interval_pos, interval_neg, raw_cutoff, raw_interval_pos, raw_interval_neg, experiment, experiment_date):
    fig = plt.figure(figsize=(10,10), dpi=100)

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.scatter(range(len(frednorm)), np.sort(frednorm), c=np.sort(d1_d2), cmap='prism_r')
    ax1.axhline(y=cutoff, color='black')
    ax1.axhline(y=interval_pos, color='green', linestyle='--')
    ax1.axhline(y=interval_neg, color='green', linestyle='--')
    ax1.set_xlabel('Cell')
    ax1.set_ylabel('Normalized Red Fluorescence')
    ax1.set_title('Norm Red Fluoresence')

    ax2.scatter(range(len(fredmean)), np.sort(fredmean), c=np.sort(d1_d2), cmap='prism_r')
    ax2.axhline(y=raw_cutoff, color='black')
    ax2.axhline(y=raw_interval_pos, color='green', linestyle='--')
    ax2.axhline(y=raw_interval_neg, color='green', linestyle='--')
    ax2.set_xlabel('Cell')
    ax2.set_ylabel('Avg Red Fluorescence')
    ax2.set_title('Avg Fluorescence')

    ax3.scatter(range(len(frednorm)), (frednorm), c=d1_d2, cmap='prism_r')
    ax3.axhline(y=cutoff, color='black')
    ax3.axhline(y=interval_pos, color='green', linestyle='--')
    ax3.axhline(y=interval_neg, color='green', linestyle='--')
    ax3.set_xlabel('Cell')
    ax3.set_ylabel('Normalized Red Fluorescence')
    ax3.set_title('Norm Red Fluoresence')

    ax4.scatter(range(len(fredmean)), (fredmean), c=d1_d2, cmap='prism_r')
    ax4.axhline(y=raw_cutoff, color='black')
    ax4.axhline(y=raw_interval_pos, color='green', linestyle='--')
    ax4.axhline(y=raw_interval_neg, color='green', linestyle='--')
    ax4.set_xlabel('Cell')
    ax4.set_ylabel('Avg Red Fluorescence')
    ax4.set_title('Avg Fluorescence')
    plt.savefig(experiment_date + '--' + experiment + '--Red_Distribution.png', bbox_inches='tight')


#%%
def _plot_found_spikes(roll_e, sav_e, d_cells, df, cells_plotted, spk_frame, Fmean_arr, Fpeak, 
                       indices, stimulations, repetition, fs, stimtime, framenum, framerate,
                       rolling, 
                       savefig = False):

    if rolling: 
        final_df = roll_e.copy()
        b3 = np.reshape(final_df, d_cells.shape)
    else: 
        final_df = sav_e.copy()
        b3 = np.reshape(sav_df, d_cells.shape)

    spike_frame = spk_frame.copy()
    Fmean = Fmean_arr.copy()

    for n in cells_plotted:
        k = len(stimulations)*repetition
        fig1 = plt.figure(figsize=(8,4), dpi = 150)
        plt.suptitle(f"ROI# {str(indices[n])}, n= {str(n)}", weight = 'bold')

        ax1 = fig1.add_subplot(211)
        ax1.plot(b3[n], color = 'k', label = 'Processed Waveform', alpha = 0.6)

        ax1.set_ylabel('F/F0', weight = 'demibold')
        ax1.set_xlabel('Frames', weight = 'demibold')
        ax1.legend()

        xcoor = framenum*np.array(range(k))
        spikeframe = spike_frame[n*k:(n+1)*k] + xcoor
        peaks = Fpeak[n*k:(n+1)*k]
        means = Fmean[n*k:(n+1)*k]

        for i in range(len(xcoor)):
            ax1.axvline(x=xcoor[i], color='green', linestyle = '--', linewidth=0.5)
            ax1.scatter(x=spikeframe[i], y=peaks[i], color='blue', marker='X', zorder = 0, alpha = 0.5)
            ax1.scatter(x=spikeframe[i], y=means[i], color='red', marker='s', zorder = 5)

        ###########################################
        # Plot out individual traces for ID'd peaks
        ###########################################

        traces = df.iloc[n*k:(n+1)*k,:]
        none = traces.loc[(traces['spike_happened']!=1) & (traces['StdDev']<0.2)]
        spikes = traces.loc[(traces['spike_happened']==1) & (traces['StdDev']<0.2)]
        noisy_none = traces.loc[(traces['spike_happened'] != 1) & (traces['StdDev']>=0.2)]
        noisy_spikes = traces.loc[(traces['spike_happened'] == 1) & (traces['StdDev']>=0.2)]

        ax2 = fig1.add_subplot(245)
        ax3 = fig1.add_subplot(246)
        ax4 = fig1.add_subplot(247)
        ax5 = fig1.add_subplot(248)

        ###
        for i in spikes.index:
            x = np.linspace(0, (len(sav_e[i])*1000/fs), framenum)
            ax2.axvline(x=stimtime-100, linestyle='--', color='red')
            ax2.axhline(y=1.2, c='k', alpha = 0.5, zorder = 0)

            spikeframe = spikes['spike_time(ms)'][i]
            peaks = spikes['dF'][i]

            ax2.scatter(spikeframe, peaks, facecolor = 'g', edgecolor = 'k')
            ax2.plot(x, sav_e[i], marker='None', linestyle='-', color = 'green', alpha = 0.3)

        ###    
        for i in noisy_spikes.index:
            x = np.linspace(0, (len(sav_e[i])*1000/fs), framenum)
            ax3.axvline(x=stimtime-100, linestyle='--', color='red')
            ax3.axhline(y=1.2, c='k', alpha = 0.5, zorder = 0)

            spikeframe = noisy_spikes['spike_time(ms)'][i]
            peaks = noisy_spikes['dF'][i]

            ax3.scatter(spikeframe, peaks, facecolor = 'r')
            ax3.plot(x, sav_e[i], marker='None', linestyle='-', color = 'red', alpha = 0.3)

        ###
        for i in none.index:
            x = np.linspace(0, len(sav_e[i])*1000/framerate, framenum)
            ax4.plot(x, sav_e[i], marker='None', linestyle='-', color = 'grey', alpha = 0.3)
            ax4.axvline(x=stimtime-100, linestyle='--', color='red')
            ax4.axhline(y=1.2, c='k', alpha = 0.5, zorder = 0)

        ###
        for i in noisy_none.index:
            x = np.linspace(0, len(sav_e[i])*1000/framerate, framenum)
            ax5.plot(x, sav_e[i], marker='None', linestyle='-', color = 'red', alpha = 0.3)
            ax5.axvline(x=stimtime-100, linestyle='--', color='red')
            ax5.axhline(y=1.2, c='k', alpha = 0.5, zorder = 0)


        ylim_a = ax1.get_ylim()[1]

        ax2.set_ylim(0, ylim_a)
        ax3.set_ylim(0, ylim_a)
        ax4.set_ylim(0, ylim_a)
        ax5.set_ylim(0, ylim_a)

        xlim_a = np.max([ax2.get_xlim()[1], ax3.get_xlim()[1], ax4.get_xlim()[1], ax5.get_xlim()[1]]) 

        ax2.set_xlim(0, xlim_a)
        ax3.set_xlim(0, xlim_a)
        ax4.set_xlim(0, xlim_a)
        ax5.set_xlim(0, xlim_a)

        ax2.set_ylabel('F/F0', weight='demibold')

        ax2.set_xlabel('Time (ms)', weight = 'demibold')
        ax3.set_xlabel('Time (ms)', weight = 'demibold')
        ax4.set_xlabel('Time (ms)', weight = 'demibold')
        ax5.set_xlabel('Time (ms)', weight = 'demibold')

        ax2.set_title('Spikes')
        ax3.set_title('Noisy Spikes', c='red')
        ax4.set_title('No Spike')
        ax5.set_title('No Spike--Noisy', c='red')
        plt.tight_layout()
        
        if savefig: 
            plt.savefig("{experiment_date} -- {experiment} -- ROI#{str(indices[n])}_sample_trace.png", bbox_inches='tight')


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


#%%


