# -*- coding: utf-8 -*-
"""
Created on Wed May 09 10:15:45 2024

@author: Carlo Vignali

@def: This is preprocessing code for fiber photometry data

@source: Modified from Photometry data preprocessing.ipynb written by Lauren Burgeno (in sourceCode folder)
"""

#%% Part 0.0 - function for processing rxPreproc output

def makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath):
    import os; import pickle
    hemis = ['left', 'right']

    if rxSite == 'left':
        hemisPreproc[rxSite] = preprocOutput  # Manipulate existing dictionary
    elif rxSite == 'right':
        hemisPreproc[rxSite] = preprocOutput  # Manipulate existing dictionary

    file_path = 'hemisPreproc_{}.pickle'.format(rxSite)
    # Writing the dictionary to a file using pickle
    with open(os.path.join(rxPath, file_path), 'wb') as file:
        pickle.dump(hemisPreproc, file)

    return hemisPreproc

#%% Part 1.0 - preprocessing section [adapted from 'photometry data preprocessing.ipynb' in sourceCode folder]
# Steps are as follows:
    # 1. Lowpass filtering to reduce noise.
    # 2. Correction for photobleaching, i.e. the slow decreace in the fluorescence signal over time. Two different methods are shown (i) subtraction of a double exponential fit (ii) highpass filtering with a very low cutoff frequency.
    # 3. Movement correction by subtracting a linear fit of the movement control channel.
    # 4. Conversion of the signal to dF/F.

def preprocessing(data405, data465, data560, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, experiment, seshNumber, thisMouse):
    
    import os
    import numpy as  np
    import matplotlib.pyplot as plt  # standard Python plotting library
    import pickle
    from scipy.signal import butter, filtfilt
    from scipy.stats import linregress
    from scipy.optimize import curve_fit
    
    os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
    import localFunctions as lf
    
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
       

    preprocOutput = {}  # Initialize an empty dictionary to store outputs
        
#%% Part 1.1 - 
    
    print('Starting preprocessing:')
    preprocOutput['isos_raw'] = data405
    preprocOutput['gACh_raw'] = data465
    preprocOutput['rdLight_raw'] = data560
    time_seconds = epoc_ticker
    preprocOutput['sampling_rate'] = fs
    preprocOutput['plotPoints'] = np.linspace(0,len(time_seconds),len(preprocOutput['isos_raw']))
    preprocOutput['beamBreak_times'] = epoc_beamBreak
    
    lowerLim_value = (min(preprocOutput['plotPoints'], key=lambda x: abs(x - 950)))
    lowerLim_index = np.argmin(np.abs(preprocOutput['plotPoints'] - 950))
    upperLim_value = (min(preprocOutput['plotPoints'], key=lambda x: abs(x - 1070)))
    upperLim_index = np.argmin(np.abs(preprocOutput['plotPoints'] - 1070))
    
    ##plot raw
    
    #RdLight1 raw
    fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_raw'], 'r', label='rdLight') #plot rdLight on left y-axis
    ax2=plt.twinx() # create a right y-axis, sharing x-axis on the same plot
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_raw'], 'b', label='isos') # plot isos on right y-axis
    
    # Plot rewards times as ticks.
    beamBreak_ticks = ax2.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), 1+np.min(preprocOutput['isos_raw'][5000:])), label='Beam Break', color='w', marker="|", mec='k')
    
    ax1.set_ylim(np.min(preprocOutput['rdLight_raw'][5000:]), np.max(preprocOutput['rdLight_raw'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][5000:]), np.max(preprocOutput['isos_raw'][5000:]))
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='r')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Raw signals')
    
    lines = plot1 + plot2 +beamBreak_ticks #line handle for legend
    labels = [l.get_label() for l in lines]  #get legend labels
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.93)) #add legend
    plt.savefig('{directory}RdLight1+Isos_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)

    #grabACh3.1+ raw
    fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_raw'], 'g', label='grabACh') # plot isos on right y-axis
    ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_raw'], 'b', label='isos') # plot isos on right y-axis
    
    # Plot rewards times as ticks.
    beamBreak_ticks = ax2.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), 1+np.min(preprocOutput['isos_raw'][5000:])), label='Beam Break', color='w', marker="|", mec='k')
    
    ax1.set_ylim(np.min(preprocOutput['gACh_raw'][5000:]), np.max(preprocOutput['gACh_raw'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][5000:]), np.max(preprocOutput['isos_raw'][5000:]))
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('grabACh Signal (V)', color='g')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Raw signals')
    
    lines = plot1 + plot2 +beamBreak_ticks #line handle for legend
    labels = [l.get_label() for l in lines]  #get legend labels
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.93)) #add legend
    plt.savefig('{directory}gACh3.1+Isos_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    
    
    #%% Part 1.2 - denoise with lowpass filter
    
    print('    Denoising with lowpass filter:')
    # Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
    b,a = butter(2, 10, btype='low', fs=preprocOutput['sampling_rate'])
    preprocOutput['rdLight_denoised'] = filtfilt(b,a, preprocOutput['rdLight_raw'])
    preprocOutput['gACh_denoised'] = filtfilt(b,a, preprocOutput['gACh_raw'])
    preprocOutput['isos_denoised'] = filtfilt(b,a, preprocOutput['isos_raw'])
    
    #RdLight1 denoised
    fig,ax1=plt.subplots()
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_denoised'], 'r', label='rdLight denoised')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], 'b', label='isos denoised')
    beamBreak_ticks = ax2.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), 1+np.min(preprocOutput['isos_denoised'][5000:])), label='Beam Break', color='w', marker="|", mec='k', ms=10)
    
    ax1.set_ylim(np.min(preprocOutput['rdLight_denoised'][5000:]), np.max(preprocOutput['rdLight_denoised'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_denoised'][5000:]), np.max(preprocOutput['isos_denoised'][5000:]))
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='r')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals')
    
    lines = plot1+plot2+beamBreak_ticks #line handle for legend
    labels = [l.get_label() for l in lines]  #get legend labels
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.92)) #add legend
    plt.savefig('{directory}denoised-rdLight1_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #grabACh denoised
    fig,ax1=plt.subplots()
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_denoised'], 'g', label='grabACh denoised')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], 'b', label='isos denoised')
    beamBreak_ticks = ax2.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), 1+np.min(preprocOutput['isos_denoised'][5000:])), label='Beam Break', color='w', marker="|", mec='k', ms=10)
    
    ax1.set_ylim(np.min(preprocOutput['gACh_denoised'][5000:]), np.max(preprocOutput['gACh_denoised'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_denoised'][5000:]), np.max(preprocOutput['isos_denoised'][5000:]))
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('grabACh Signal (V)', color='g')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals')
    
    lines = plot1+plot2+beamBreak_ticks #line handle for legend
    labels = [l.get_label() for l in lines]  #get legend labels
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.98, 0.92)) #add legend
    plt.savefig('{directory}denoised-gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)   
    
    #%% Part 1.3 - plot smaller window of post-denoise data
    
    #RdLight1 denoised (small window)
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_raw'], color='r', alpha=0.3, label='rdLight raw')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_raw'], color='b', alpha=0.3, label='isos raw') 
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_denoised'], color='r', label='rdLight denoised') 
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], color='b', label='isos denoised') 
    beamBreak_ticks = ax1.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['rdLight_denoised'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='r')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals')
    
    lines = plot1+plot2 + plot3 + plot4 + beamBreak_ticks
    labels = [l.get_label() for l in lines]
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
    ax1.set_xlim(lowerLim_value, upperLim_value) # 60 sec window
    ax1.set_ylim(-5+np.min(preprocOutput['rdLight_raw'][lowerLim_index:upperLim_index]), np.max(preprocOutput['rdLight_raw'][lowerLim_index:upperLim_index]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][lowerLim_index:upperLim_index]), 5+np.max(preprocOutput['isos_raw'][lowerLim_index:upperLim_index]))
    plt.savefig('{directory}denoised60sWindow-rdLight1_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #grabACh denoised (small window)
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_raw'], color='g', alpha=0.3, label='grabACh raw')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_raw'], color='b', alpha=0.3, label='isos raw') 
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_denoised'], color='g', label='grabACh denoised') 
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], color='b', label='isos denoised') 
    beamBreak_ticks = ax1.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['gACh_denoised'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('grabACh Signal (V)', color='g')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals')
    
    lines = plot1+plot2 + plot3 + plot4 + beamBreak_ticks
    labels = [l.get_label() for l in lines]
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
    ax1.set_xlim(lowerLim_value, upperLim_value) # 60 sec window
    ax1.set_ylim(-5+np.min(preprocOutput['gACh_raw'][lowerLim_index:upperLim_index]), np.max(preprocOutput['gACh_raw'][lowerLim_index:upperLim_index]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][lowerLim_index:upperLim_index]), 5+np.max(preprocOutput['isos_raw'][lowerLim_index:upperLim_index]))
    plt.savefig('{directory}denoised60sWindow-gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #denoised comp
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_raw'], color='g', alpha=0.3, label='grabACh raw')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_raw'], color='r', alpha=0.3, label='rdLight raw') 
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_denoised'], color='g', label='grabACh denoised') 
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_denoised'], color='r', label='rdLight denoised') 
    beamBreak_ticks = ax1.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['gACh_denoised'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('grabACh Signal (V)', color='g')
    ax2.set_ylabel('rdLight Signal (V)', color='r')
    ax1.set_title('Denoised signals')
    
    lines = plot1+plot2 + plot3 + plot4 + beamBreak_ticks
    labels = [l.get_label() for l in lines]
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.93, 0.99))
    ax1.set_xlim(lowerLim_value, upperLim_value) # 60 sec window
    ax1.set_ylim(-5+np.min(preprocOutput['gACh_raw'][lowerLim_index:upperLim_index]), np.max(preprocOutput['gACh_raw'][lowerLim_index:upperLim_index]))
    ax2.set_ylim(np.min(preprocOutput['rdLight_raw'][lowerLim_index:upperLim_index]), 5+np.max(preprocOutput['rdLight_raw'][lowerLim_index:upperLim_index]))
    plt.savefig('{directory}denoised60sWindow-rdLight+gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #%% Part 1.4 - photobleach correction option #1 > double exponential
    
    print('    Correcting photobleaching with double exponential:')
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
    
    def weighted_tau(amp_fast, amp_slow, tau_slow, tau_multiplier):
        tau_fast = tau_slow * tau_multiplier
        return (amp_fast * tau_fast + amp_slow * tau_slow) / (amp_fast + amp_slow)
    def initial_slope(const, amp_fast, amp_slow, tau_slow, tau_multiplier):
        tau_fast = tau_slow * tau_multiplier
        return - (amp_slow / tau_slow + amp_fast / tau_fast)
    ## !! correction for rdLight !!
    
    # Fit curve to rdLight signal.
    max_sig = np.max(preprocOutput['rdLight_denoised'])
    inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
    bounds = ([0      , 0      , 0      , 600  , 0],
              [max_sig, max_sig, max_sig, 36000, 1])
    preprocOutput['rdLight_parms'], parm_cov = curve_fit(double_exponential, preprocOutput['plotPoints'], preprocOutput['rdLight_denoised'], 
                                      p0=inital_params, bounds=bounds, maxfev=1000)
    preprocOutput['rdLight_expfit'] = double_exponential(preprocOutput['plotPoints'], *preprocOutput['rdLight_parms'])
    preprocOutput['rdLight_effTau'] = weighted_tau(*preprocOutput['rdLight_parms'][1:])
    preprocOutput['rdLight_decaySlope'] = initial_slope(*preprocOutput['rdLight_parms'])
    
    # Fit curve to isos signal.
    max_sig = np.max(preprocOutput['isos_denoised'])
    inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
    bounds = ([0      , 0      , 0      , 600  , 0],
              [max_sig, max_sig, max_sig, 36000, 1])
    preprocOutput['isos_parms'], parm_cov = curve_fit(double_exponential, preprocOutput['plotPoints'], preprocOutput['isos_denoised'], 
                                      p0=inital_params, bounds=bounds, maxfev=1000)
    preprocOutput['isos_expfit'] = double_exponential(preprocOutput['plotPoints'], *preprocOutput['isos_parms'])
    preprocOutput['isos_effTau'] = weighted_tau(*preprocOutput['isos_parms'][1:])
    preprocOutput['isos_decaySlope'] = initial_slope(*preprocOutput['isos_parms'])
    
    #plot fits over denoised data
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_denoised'], 'r', label='rdLight')
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_expfit'], color='k', linewidth=1.5, label='Exponential fit') 
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], color='b', label='isos') 
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_expfit'],color='k', linewidth=1.5) 
    
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='r')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals with double exponential fits')
    
    lines = plot1 + plot2 + plot3
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right'); 
    ax1.set_ylim(np.min(preprocOutput['rdLight_raw'][5000:]), np.max(preprocOutput['rdLight_raw'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][5000:]), np.max(preprocOutput['isos_raw'][5000:]))
    plt.savefig('{directory}bleachExponentialFit-rdLight1_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    ##subtract exponentials fit to data
    
    preprocOutput['rdLight_detrended'] = preprocOutput['rdLight_denoised'] - preprocOutput['rdLight_expfit']
    preprocOutput['isos_detrended'] = preprocOutput['isos_denoised'] - preprocOutput['isos_expfit']
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_detrended'], 'r', label='rdLight')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_detrended'], color='b', label='isos') 
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='r')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Bleaching Correction by Double Exponential Fit')
    
    lines = plot1+plot2 
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right'); 
    ax1.set_xlim(lowerLim_value, upperLim_value) # 60 sec window
    ax1.set_ylim(-3+np.min(preprocOutput['rdLight_detrended'][5000:]), np.max(preprocOutput['rdLight_detrended'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_detrended'][5000:]), 3+np.max(preprocOutput['isos_detrended'][5000:]))
    plt.savefig('{directory}bleachCorrection-rdLight1_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    
    ## !! correction for gACh !!
    
    # Fit curve to gACh signal.
    max_sig = np.max(preprocOutput['gACh_denoised'])
    inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
    bounds = ([0      , 0      , 0      , 600  , 0],
              [max_sig, max_sig, max_sig, 36000, 1])
    preprocOutput['gACh_parms'], parm_cov = curve_fit(double_exponential, preprocOutput['plotPoints'], preprocOutput['gACh_denoised'], 
                                      p0=inital_params, bounds=bounds, maxfev=1000)
    preprocOutput['gACh_expfit'] = double_exponential(preprocOutput['plotPoints'], *preprocOutput['gACh_parms'])
    preprocOutput['gACh_effTau'] = weighted_tau(*preprocOutput['gACh_parms'][1:])
    preprocOutput['gACh_decaySlope'] = initial_slope(*preprocOutput['gACh_parms'])
    
    # Fit curve to isos signal.
        #calculations are above
    
    #plot fits over denoised data
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_denoised'], 'g', label='gACh')
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_expfit'], color='k', linewidth=1.5, label='Exponential fit') 
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_denoised'], color='b', label='isos') 
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_expfit'],color='k', linewidth=1.5) 
    
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('gACh Signal (V)', color='g')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Denoised signals with double exponential fits')
    
    lines = plot1 + plot2 + plot3
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right'); 
    ax1.set_ylim(np.min(preprocOutput['gACh_raw'][5000:]), np.max(preprocOutput['gACh_raw'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_raw'][5000:]), np.max(preprocOutput['isos_raw'][5000:]))
    plt.savefig('{directory}bleachExponentialFit-gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    ##subtract exponentials fit to data
    
    preprocOutput['gACh_detrended'] = preprocOutput['gACh_denoised'] - preprocOutput['gACh_expfit']
    preprocOutput['isos_detrended'] = preprocOutput['isos_denoised'] - preprocOutput['isos_expfit']
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_detrended'], 'g', label='gACh')
    ax2=plt.twinx()
    plot2=ax2.plot(preprocOutput['plotPoints'], preprocOutput['isos_detrended'], color='b', label='isos') 
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('gACh Signal (V)', color='g')
    ax2.set_ylabel('isos Signal (V)', color='b')
    ax1.set_title('Bleaching Correction by Double Exponential Fit')
    
    lines = plot1+plot2 
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right'); 
    #ax1.set_xlim(lowerLim_value, upperLim_value) # 60 sec window
    ax1.set_ylim(-5+np.min(preprocOutput['gACh_detrended'][5000:]), np.max(preprocOutput['gACh_detrended'][5000:]))
    ax2.set_ylim(np.min(preprocOutput['isos_detrended'][5000:]), 5+np.max(preprocOutput['isos_detrended'][5000:]))
    plt.savefig('{directory}bleachCorrection-gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #plt.close('all')
    #%% Part 1.5 - motion correction
    
    print('    Motion correction:')
    ## !! for RdLight1 signal !!
    
    #regression between isos and fluorophore
    slope, intercept, r_value, p_value, std_err = linregress(x=preprocOutput['isos_detrended'], y=preprocOutput['rdLight_detrended'])
    
    fig,ax1=plt.subplots() 
    plt.scatter(preprocOutput['isos_detrended'][::5], preprocOutput['rdLight_detrended'][::5],alpha=0.1, marker='.')
    x = np.array(plt.xlim())
    plt.plot(x, intercept+slope*x)
    plt.xlabel('isos')
    plt.ylabel('rdLight')
    plt.title('isos - rdLight correlation.')
    
    print('Slope    : {:.3f}'.format(slope))
    print('R-squared: {:.3f}'.format(r_value**2))
    
    ##plot
    
    #estimated motion-related signal in fluorophore
    preprocOutput['rdLight_est_motion'] = intercept + slope * preprocOutput['isos_detrended']
    #subtraction of motion-related signal
    preprocOutput['rdLight_corrected'] = preprocOutput['rdLight_detrended'] - preprocOutput['rdLight_est_motion']
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_detrended'], 'k' , label='rdLight - pre motion correction', alpha=0.5)
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_corrected'], 'r', label='rdLight - motion corrected', alpha=0.5)
    ax2=plt.twinx()
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['rdLight_est_motion'] - 0.05, 'y', label='estimated motion')
    beamBreak_ticks = ax2.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['rdLight_est_motion'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('rdLight Signal (V)', color='g')
    ax1.set_title('Motion Correction')
    
    lines = plot1+plot3+plot4 + beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    ax1.set_xlim(lowerLim_value, upperLim_value)  # 60 sec window
    ax1.set_ylim(-4+np.min(preprocOutput['rdLight_detrended'][lowerLim_index:upperLim_index]), 3+np.max(preprocOutput['rdLight_detrended'][lowerLim_index:upperLim_index]))
    ax2.set_ylim(-1+np.min(preprocOutput['rdLight_est_motion'][lowerLim_index:upperLim_index]), 6+np.max(preprocOutput['rdLight_est_motion'][lowerLim_index:upperLim_index]))
    plt.savefig('{directory}motionCorrected-rdLight1_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    
    ## !! for grabACh signal !!
    
    slope, intercept, r_value, p_value, std_err = linregress(x=preprocOutput['isos_detrended'], y=preprocOutput['gACh_detrended'])
    
    fig,ax1=plt.subplots() 
    plt.scatter(preprocOutput['isos_detrended'][::5], preprocOutput['gACh_detrended'][::5],alpha=0.1, marker='.')
    x = np.array(plt.xlim())
    plt.plot(x, intercept+slope*x)
    plt.xlabel('isos')
    plt.ylabel('gACh')
    plt.title('isos - gACh correlation.')
    
    print('Slope    : {:.3f}'.format(slope))
    print('R-squared: {:.3f}'.format(r_value**2))
    
    ##plot
    
    preprocOutput['gACh_est_motion'] = intercept + slope * preprocOutput['isos_detrended']
    preprocOutput['gACh_corrected'] = preprocOutput['gACh_detrended'] - preprocOutput['gACh_est_motion']
    
    fig,ax1=plt.subplots()  
    plot1=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_detrended'], 'k' , label='gACh - pre motion correction', alpha=0.5)
    plot3=ax1.plot(preprocOutput['plotPoints'], preprocOutput['gACh_corrected'], 'g', label='gACh - motion corrected', alpha=0.5)
    ax2=plt.twinx()
    plot4=ax2.plot(preprocOutput['plotPoints'], preprocOutput['gACh_est_motion'] - 0.05, 'y', label='estimated motion')
    beamBreak_ticks = ax1.plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['gACh_est_motion'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('gACh Signal (V)', color='g')
    ax1.set_title('Motion Correction')
    
    lines = plot1+plot3+plot4 + beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    ax1.set_xlim(lowerLim_value, upperLim_value)  # 60 sec window
    ax1.set_ylim(-4+np.min(preprocOutput['gACh_corrected'][lowerLim_index:upperLim_index]), 3+np.max(preprocOutput['gACh_corrected'][lowerLim_index:upperLim_index]))
    ax2.set_ylim(-1+np.min(preprocOutput['gACh_est_motion'][lowerLim_index:upperLim_index]), 6+np.max(preprocOutput['gACh_est_motion'][lowerLim_index:upperLim_index]))
    plt.savefig('{directory}motionCorrected-gACh_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #%% Part 1.6 - normalize dF/F
    
    print('    Plotting DFoF:')
    preprocOutput['rdLight_dF_F'] = 100*preprocOutput['rdLight_corrected']/preprocOutput['rdLight_expfit']
    preprocOutput['gACh_dF_F'] = 100*preprocOutput['gACh_corrected']/preprocOutput['gACh_expfit']
    
    fig,axs=plt.subplots(2,1)  
    plot1=axs[0].plot(preprocOutput['plotPoints'], preprocOutput['rdLight_dF_F'], 'r', label='rdLight dF/F')
    beamBreak_ticks = axs[0].plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['rdLight_dF_F'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    axs[0].set_xlabel('Time (seconds)')
    axs[0].set_ylabel('dF/F (%)')
    axs[0].set_title('dF/F for rDLight')
    
    lines = plot1+ beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = axs[0].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    axs[0].set_xlim(lowerLim_value, upperLim_value)
    axs[0].set_ylim(-5,5);
    
    plot2=axs[1].plot(preprocOutput['plotPoints'], preprocOutput['gACh_dF_F'], 'g', label='gACh dF/F')
    beamBreak_ticks = axs[1].plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['gACh_dF_F'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    axs[1].set_xlabel('Time (seconds)')
    axs[1].set_ylabel('dF/F (%)')
    axs[1].set_title('dF/F for gACh')
    
    lines = plot2+ beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = axs[1].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    axs[1].set_xlim(lowerLim_value, upperLim_value)
    axs[1].set_ylim(-5,5);
    plt.savefig('{directory}normalizedDFoF_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    # !! for isos channel !!
    preprocOutput['isos_dF_F'] = 100*preprocOutput['isos_detrended']/preprocOutput['isos_expfit']
    
    fig,axs=plt.subplots(2,1)  
    plot1=axs[0].plot(preprocOutput['plotPoints'], preprocOutput['isos_dF_F'], 'b', label='isos dF/F')
    beamBreak_ticks = axs[0].plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['isos_dF_F'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    axs[0].set_xlabel('Time (seconds)')
    axs[0].set_ylabel('dF/F (%)')
    axs[0].set_title('dF/F for isos')
    
    lines = plot1+ beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = axs[0].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    axs[0].set_xlim(lowerLim_value, upperLim_value)
    axs[0].set_ylim(-5,5);
    
    plot2=axs[1].plot(preprocOutput['plotPoints'], preprocOutput['gACh_dF_F'], 'g', label='gACh dF/F')
    beamBreak_ticks = axs[1].plot(preprocOutput['beamBreak_times'], np.full(np.size(preprocOutput['beamBreak_times']), np.min(preprocOutput['gACh_dF_F'][lowerLim_index:upperLim_index])-1), label='Beam Break',color='w', marker="|", mfc='k', mec='k', ms=8)
    
    axs[1].set_xlabel('Time (seconds)')
    axs[1].set_ylabel('dF/F (%)')
    axs[1].set_title('dF/F for gACh')
    
    lines = plot2+ beamBreak_ticks
    labels = [l.get_label() for l in lines]  
    legend = axs[1].legend(lines, labels, loc='upper right', bbox_to_anchor=(0.95, 0.98))
    
    axs[1].set_xlim(lowerLim_value, upperLim_value)
    axs[1].set_ylim(-5,5);
    plt.savefig('{directory}normalizedDFoF_ofIsos_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=rxSite, exp=experiment, sesh=seshNumber,
mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    #%% Part 2.0 - return / save
    
    print(' Done.')
    # File path
    file_path = 'preprocOutput_{}.pickle'.format(rxSite)
    # Writing the dictionary to a file using pickle
    with open(os.path.join(rxPath, file_path), 'wb') as file:
        pickle.dump(preprocOutput, file)
    
    return preprocOutput