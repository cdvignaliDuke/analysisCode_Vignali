# -*- coding: utf-8 -*-
"""
Created on Wed May 09 10:15:45 2024

@author: Carlo Vignali

@def: This is alignment code for fiber photometry data (aligning data to set epocs)

@source: 
"""

def alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex):

#%% Part 0.0 - preprocessing section [adapted from 'photometry data preprocessing.ipynb' in sourceCode folder]
# Steps are as follows:
    # 1. Lowpass filtering to reduce noise.
    # 2. Correction for photobleaching, i.e. the slow decreace in the fluorescence signal over time. Two different methods are shown (i) subtraction of a double exponential fit (ii) highpass filtering with a very low cutoff frequency.
    # 3. Movement correction by subtracting a linear fit of the movement control channel.
    # 4. Conversion of the signal to dF/F.
    
    import os
    import numpy as  np
    import matplotlib.pyplot as plt  # standard Python plotting library
    import pickle
    from matplotlib.lines import Line2D
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
             
    #%% Part 1.0 - format data to align to epocs
    print('Aligning rx data to bx events:')
    
    alignmentOutput={}; alignedData={}; 
    
    if rxSite == 'left':
        sites=['left']
        alignmentOutput['ACh_dFoF']={}; alignmentOutput['DA_dFoF']={}; alignmentOutput['isos_dFoF']={};
        for s in range(len(sites)):
            alignmentOutput['ACh_dFoF'][s] = hemisPreproc[sites[s]]['gACh_dF_F']
            alignmentOutput['DA_dFoF'][s] = hemisPreproc[sites[s]]['rdLight_dF_F']
            alignmentOutput['isos_dFoF'][s] = hemisPreproc[sites[s]]['isos_dF_F']
            alignmentOutput['plotPoints'] = hemisPreproc[sites[s]]['plotPoints']
    elif rxSite == 'right':
        sites=['right']
        alignmentOutput['ACh_dFoF']={}; alignmentOutput['DA_dFoF']={}; alignmentOutput['isos_dFoF']={};
        for s in range(len(sites)):
            alignmentOutput['ACh_dFoF'][s] = hemisPreproc[sites[s]]['gACh_dF_F']
            alignmentOutput['DA_dFoF'][s] = hemisPreproc[sites[s]]['rdLight_dF_F']
            alignmentOutput['isos_dFoF'][s] = hemisPreproc[sites[s]]['isos_dF_F']
            alignmentOutput['plotPoints'] = hemisPreproc[sites[s]]['plotPoints']
    elif rxSite == 'both':
        sites=['left','right']
        alignmentOutput['ACh_dFoF']={}; alignmentOutput['DA_dFoF']={}; alignmentOutput['isos_dFoF']={};
        for s in range(len(sites)):
            alignmentOutput['ACh_dFoF'][s] = hemisPreproc[sites[s]]['gACh_dF_F']
            alignmentOutput['DA_dFoF'][s] = hemisPreproc[sites[s]]['rdLight_dF_F']
            alignmentOutput['isos_dFoF'][s] = hemisPreproc[sites[s]]['isos_dF_F']
            alignmentOutput['plotPoints'] = hemisPreproc['left']['plotPoints']
    
            
    if len(epoc_beamBreak) <= 1: 
        print('    Insufficient beamBreak for alignment:')
        
        alignmentOutput['index_trialWithStartCue'] = {}; alignmentOutput['index_trialWithTone'] = {}; alignmentOutput['value_trialWithStartCue'] = {}; alignmentOutput['value_trialWithTone'] = {}; 
        alignmentOutput['nTrials'] = [val for val in range(len(epoc_onsetITI)) if epoc_onsetITI[val] < 3600];
        
        for trial in alignmentOutput['nTrials']:
            if epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] if np.any(epoc_onsetCues > epoc_onsetITI[trial]) else False:
                alignmentOutput['index_trialWithStartCue'][s][trial] = True
                alignmentOutput['value_trialWithStartCue'][s][trial] = epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] 
            if np.where(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5)[0][-1] if np.any(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5) else False:
                alignmentOutput['index_trialWithTone'][s][trial] = True
                alignmentOutput['value_trialWithTone'][s][trial] = epoc_onsetCues[np.where(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5)[0][-1]]
        
    elif len(epoc_beamBreak) > 1:
        print('    Align to beam break:')
    
        alignmentOutput['index_trialWithTone'] = {}; alignmentOutput['index_trialWithStartCue'] = {}; alignmentOutput['index_trialWithResponse'] = {}; alignmentOutput['index_trialWithTooEarlyResponse'] = {}; alignmentOutput['value_trialWithResponse'] = {}; alignmentOutput['value_trialWithTooEarlyResponse'] = {}; alignmentOutput['value_trialWithTone'] = {}; alignmentOutput['value_trialWithStartCue'] = {}
        for s in range(len(sites)):
            alignmentOutput['index_trialWithStartCue'][s] = np.zeros(len(epoc_onsetITI),dtype=bool); 
            alignmentOutput['index_trialWithTone'][s] = np.zeros(len(epoc_onsetITI),dtype=bool); 
            alignmentOutput['index_trialWithResponse'][s] = np.zeros(len(epoc_onsetITI),dtype=bool); 
            alignmentOutput['index_trialWithTooEarlyResponse'][s] = np.zeros(len(epoc_onsetITI),dtype=bool)
            alignmentOutput['value_trialWithStartCue'][s] = np.zeros(len(epoc_onsetITI)); 
            alignmentOutput['value_trialWithTone'][s] = np.zeros(len(epoc_onsetITI)); 
            alignmentOutput['value_trialWithResponse'][s] = np.zeros(len(epoc_onsetITI)); 
            alignmentOutput['value_trialWithTooEarlyResponse'][s] = np.zeros(len(epoc_onsetITI))
            alignmentOutput['nTrials'] = [val for val in range(len(epoc_onsetITI)) if epoc_onsetITI[val] < 3600];
                
            for trial in alignmentOutput['nTrials']:
                if keyDataSheet['exp_seshPhase'] > 1:
                    if epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] if np.any(epoc_onsetCues > epoc_onsetITI[trial]) else False:
                        alignmentOutput['index_trialWithStartCue'][s][trial] = True
                        alignmentOutput['value_trialWithStartCue'][s][trial] = epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] 
                    if not np.isnan(keyDataSheet['accOfResponse']).any():
                        if trial < len(keyDataSheet['accOfResponse']):
                            if keyDataSheet['accOfResponse'][trial] == 5:
                                alignmentOutput['index_trialWithTone'][s][trial] = False
                                alignmentOutput['value_trialWithTone'][s][trial] = 0
                            else: 
                                if np.where(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5)[0][-1] if np.any(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5) else False:
                                    alignmentOutput['index_trialWithTone'][s][trial] = True
                                    alignmentOutput['value_trialWithTone'][s][trial] = epoc_onsetCues[np.where(np.abs(epoc_onsetCues - alignmentOutput['value_trialWithStartCue'][s][trial]) <= 5)[0][-1]]
                    if alignmentOutput['value_trialWithTone'][s][trial] != 0:
                        if epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithTone'][s][trial])]:  
                            potentialBeamBreak = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithTone'][s][trial])]
                            if np.any(potentialBeamBreak > alignmentOutput['value_trialWithTone'][s][trial]) and np.any(potentialBeamBreak < alignmentOutput['value_trialWithTone'][s][trial]+5):
                                alignmentOutput['index_trialWithResponse'][s][trial] = True
                                alignmentOutput['value_trialWithResponse'][s][trial] = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithTone'][s][trial])]
                    if epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]:  
                        potentialBeamBreak = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]
                        if np.any(potentialBeamBreak > alignmentOutput['value_trialWithStartCue'][s][trial]) and np.any(potentialBeamBreak < alignmentOutput['value_trialWithStartCue'][s][trial]+(alignmentOutput['value_trialWithTone'][s][trial]-alignmentOutput['value_trialWithStartCue'][s][trial])):
                            alignmentOutput['index_trialWithTooEarlyResponse'][s][trial] = True
                            alignmentOutput['value_trialWithTooEarlyResponse'][s][trial] = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]
                elif keyDataSheet['exp_seshPhase'] == 1:
                    if epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] if np.any(epoc_onsetCues > epoc_onsetITI[trial]) else False:
                        alignmentOutput['index_trialWithStartCue'][s][trial] = True
                        alignmentOutput['value_trialWithStartCue'][s][trial] = epoc_onsetCues[np.argmax(epoc_onsetCues > epoc_onsetITI[trial])] 
                    alignmentOutput['value_trialWithTone'][s][trial] = 0
                    if epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]:  
                        potentialBeamBreak = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]
                        if np.any(potentialBeamBreak > alignmentOutput['value_trialWithStartCue'][s][trial]): #and np.any(potentialBeamBreak < alignmentOutput['value_trialWithStartCue'][s][trial]+5):
                            alignmentOutput['index_trialWithResponse'][s][trial] = True
                            alignmentOutput['value_trialWithResponse'][s][trial] = epoc_beamBreak[np.argmax(epoc_beamBreak > alignmentOutput['value_trialWithStartCue'][s][trial])]
                    
    alignmentOutput['preEpocWindow'] = 6 #seconds
    alignmentOutput['postEpocWindow'] = 11 
    alignedData['cueAlign_ACh_dFoF']={}; alignedData['cueAlign_DA_dFoF']={}; alignedData['cueAlign_isos_dFoF']={}; 
    alignedData['beamAlign_ACh_dFoF']={}; alignedData['beamAlign_DA_dFoF']={}; alignedData['beamAlign_isos_dFoF']={}; 
    for s in range(len(sites)):
        alignedData['cueAlign_ACh_dFoF'][s]=[]; alignedData['cueAlign_DA_dFoF'][s]=[]; alignedData['cueAlign_isos_dFoF'][s]=[]; 
        alignedData['beamAlign_ACh_dFoF'][s]=[]; alignedData['beamAlign_DA_dFoF'][s]=[]; alignedData['beamAlign_isos_dFoF'][s]=[]; 
        for trial in range(len(alignmentOutput['value_trialWithStartCue'][s])):
            #find index
            eventIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-alignmentOutput['value_trialWithStartCue'][s][trial]))
            preEpocIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-(alignmentOutput['value_trialWithStartCue'][s][trial]-alignmentOutput['preEpocWindow'])))
            postEpocIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-(alignmentOutput['value_trialWithStartCue'][s][trial]+alignmentOutput['postEpocWindow'])))
            #align dFoF
            alignedData['cueAlign_isos_dFoF'][s].append(alignmentOutput['isos_dFoF'][s][preEpocIndex:postEpocIndex])
            alignedData['cueAlign_ACh_dFoF'][s].append(alignmentOutput['ACh_dFoF'][s][preEpocIndex:postEpocIndex])
            alignedData['cueAlign_DA_dFoF'][s].append(alignmentOutput['DA_dFoF'][s][preEpocIndex:postEpocIndex])
                
        for trial in range(len(alignmentOutput['value_trialWithResponse'][s])):
            #find index
            eventIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-alignmentOutput['value_trialWithResponse'][s][trial]))
            preEpocIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-(alignmentOutput['value_trialWithResponse'][s][trial]-alignmentOutput['preEpocWindow'])))
            postEpocIndex = np.argmin(np.abs(alignmentOutput['plotPoints']-(alignmentOutput['value_trialWithResponse'][s][trial]+alignmentOutput['postEpocWindow'])))
            #align dFoF
            alignedData['beamAlign_isos_dFoF'][s].append(alignmentOutput['isos_dFoF'][s][preEpocIndex:postEpocIndex])
            alignedData['beamAlign_ACh_dFoF'][s].append(alignmentOutput['ACh_dFoF'][s][preEpocIndex:postEpocIndex])
            alignedData['beamAlign_DA_dFoF'][s].append(alignmentOutput['DA_dFoF'][s][preEpocIndex:postEpocIndex])     
        
    #if rxSite == 'both':
    alignmentOutput['trialOutcome']={}; alignmentOutput['trialType']={}; alignmentOutput['trialType_calc']={}; alignmentOutput['beamAlign_trialOutcome']={}; alignmentOutput['beamAlign_trialType']={}; alignmentOutput['beamAlign_trialType_calc']={};
    for s in range(len(sites)):
        for i in range(len(alignedData['cueAlign_DA_dFoF'])):
            #session ends before last trial terminates    
            alignedData['cueAlign_ACh_dFoF'][s].pop(); alignedData['cueAlign_DA_dFoF'][s].pop(); alignedData['beamAlign_ACh_dFoF'][s].pop(); alignedData['beamAlign_DA_dFoF'][s].pop(); alignedData['cueAlign_isos_dFoF'][s].pop(); alignedData['beamAlign_isos_dFoF'][s].pop();          
            #
            alignmentOutput['minFrames'] = min(arr.shape[0] for arr in alignedData['cueAlign_ACh_dFoF'][s])
            #
            alignedData['cueAlign_isos_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['cueAlign_isos_dFoF'][s]]
            alignedData['cueAlign_ACh_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['cueAlign_ACh_dFoF'][s]]
            alignedData['cueAlign_DA_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['cueAlign_DA_dFoF'][s]]
            #
            alignedData['beamAlign_isos_dFoF'][s] = [val for val, keep in zip(alignedData['beamAlign_isos_dFoF'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
            alignedData['beamAlign_ACh_dFoF'][s] = [val for val, keep in zip(alignedData['beamAlign_ACh_dFoF'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
            alignedData['beamAlign_DA_dFoF'][s] = [val for val, keep in zip(alignedData['beamAlign_DA_dFoF'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
            #
            alignedData['beamAlign_isos_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['beamAlign_isos_dFoF'][s]]
            alignedData['beamAlign_ACh_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['beamAlign_ACh_dFoF'][s]]
            alignedData['beamAlign_DA_dFoF'][s] = [arr[:alignmentOutput['minFrames']] for arr in alignedData['beamAlign_DA_dFoF'][s]]
            if not np.isnan(keyDataSheet['accOfResponse']).any():
                alignmentOutput['trialOutcome'][s] = keyDataSheet['accOfResponse'][0:len(alignedData['cueAlign_ACh_dFoF'][s])] 
                alignmentOutput['trialType'][s] = keyDataSheet['trialType'][0:len(alignedData['cueAlign_ACh_dFoF'][s])] 
                alignmentOutput['trialType_calc'][s] = keyDataSheet['trialType_calc'][0:len(alignedData['cueAlign_ACh_dFoF'][s])] 
                alignmentOutput['beamAlign_trialOutcome'][s] = [val for val, keep in zip(alignmentOutput['trialOutcome'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
                alignmentOutput['beamAlign_trialType'][s] = [val for val, keep in zip(alignmentOutput['trialType'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
                alignmentOutput['beamAlign_trialType_calc'][s] = [val for val, keep in zip(alignmentOutput['trialType_calc'][s], alignmentOutput['index_trialWithResponse'][s][:-1]) if keep]
    
    if not np.isnan(keyDataSheet['accOfResponse']).any():
        for s in range(len(sites)):
            #ensure index lengths match data (on off chance medPC script collected a trial at the very end of the session - near 3600 seconds)    
            if len(alignedData['cueAlign_ACh_dFoF'][s]) < len(alignmentOutput['index_trialWithStartCue'][s]):
                alignmentOutput['index_trialWithResponse'][s] = alignmentOutput['index_trialWithResponse'][s][:-1]
                alignmentOutput['index_trialWithStartCue'][s] = alignmentOutput['index_trialWithStartCue'][s][:-1]
                alignmentOutput['index_trialWithTone'][s] = alignmentOutput['index_trialWithTone'][s][:-1]
                alignmentOutput['index_trialWithTooEarlyResponse'][s] = alignmentOutput['index_trialWithTooEarlyResponse'][s][:-1]
               
            if len(alignmentOutput['trialOutcome'][s]) < len(alignedData['cueAlign_ACh_dFoF'][s]):
                alignedData['cueAlign_isos_dFoF'][s] = alignedData['cueAlign_isos_dFoF'][s][:-(len(alignedData['cueAlign_isos_dFoF'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignedData['cueAlign_ACh_dFoF'][s] = alignedData['cueAlign_ACh_dFoF'][s][:-(len(alignedData['cueAlign_ACh_dFoF'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignedData['cueAlign_DA_dFoF'][s] = alignedData['cueAlign_DA_dFoF'][s][:-(len(alignedData['cueAlign_DA_dFoF'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignmentOutput['index_trialWithResponse'][s] = alignmentOutput['index_trialWithResponse'][s][:-(len(alignmentOutput['index_trialWithResponse'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignmentOutput['index_trialWithStartCue'][s] = alignmentOutput['index_trialWithStartCue'][s][:-(len(alignmentOutput['index_trialWithStartCue'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignmentOutput['index_trialWithTone'][s] =alignmentOutput['index_trialWithTone'][s][:-(len(alignmentOutput['index_trialWithTone'][s])-len(alignmentOutput['trialOutcome'][s]))]
                alignmentOutput['index_trialWithTooEarlyResponse'][s] =alignmentOutput['index_trialWithTooEarlyResponse'][s][:-(len(alignmentOutput['index_trialWithTooEarlyResponse'][s])-len(alignmentOutput['trialOutcome'][s]))]
       
    #%% Part 1.1.0 - mean + std and separation based on trial outcome    

    print('    Creating session averaged data arrays')
   
    alignedData['cueAlign_ACh_stack'] = {}; alignedData['cueAlign_DA_stack'] = {}; alignedData['cueAlign_isos_stack'] = {}; alignedData['beamAlign_DA_stack'] = {}; alignedData['beamAlign_ACh_stack'] = {}; alignedData['beamAlign_isos_stack'] = {};
    for s in range(len(sites)):
        alignmentOutput['xAxisPoints'] = np.linspace(-alignmentOutput['preEpocWindow'],alignmentOutput['postEpocWindow'],alignmentOutput['minFrames'])#np.argmin(np.abs(alignmentOutput['plotPoints']-(0+alignmentOutput['preEpocWindow'])))+np.argmin(np.abs(alignmentOutput['plotPoints']-(0+alignmentOutput['postEpocWindow']))))
        
        alignedData['cueAlign_isos_dFoF'][s] = [alignedData['cueAlign_isos_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['cueAlign_isos_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['cueAlign_isos_dFoF'][s][t] for t in range(len(alignedData['cueAlign_isos_dFoF'][s]))]
        alignedData['cueAlign_ACh_dFoF'][s] = [alignedData['cueAlign_ACh_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['cueAlign_ACh_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['cueAlign_ACh_dFoF'][s][t] for t in range(len(alignedData['cueAlign_ACh_dFoF'][s]))]
        alignedData['cueAlign_DA_dFoF'][s] = [alignedData['cueAlign_DA_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['cueAlign_DA_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['cueAlign_DA_dFoF'][s][t] for t in range(len(alignedData['cueAlign_DA_dFoF'][s]))]
        alignedData['beamAlign_isos_dFoF'][s] = [alignedData['beamAlign_isos_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['beamAlign_isos_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['beamAlign_isos_dFoF'][s][t] for t in range(len(alignedData['beamAlign_isos_dFoF'][s]))]
        alignedData['beamAlign_ACh_dFoF'][s] = [alignedData['beamAlign_ACh_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['beamAlign_ACh_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['beamAlign_ACh_dFoF'][s][t] for t in range(len(alignedData['beamAlign_ACh_dFoF'][s]))]
        alignedData['beamAlign_DA_dFoF'][s] = [alignedData['beamAlign_DA_dFoF'][s][t][:alignmentOutput['minFrames']] if len(alignedData['beamAlign_DA_dFoF'][s][t]) > alignmentOutput['minFrames'] else alignedData['beamAlign_DA_dFoF'][s][t] for t in range(len(alignedData['beamAlign_DA_dFoF'][s]))]
        
        fullTrial_cue=[t for t, arr in enumerate(alignedData['cueAlign_isos_dFoF'][s]) if len(arr) == max(len(a) for a in alignedData['cueAlign_isos_dFoF'][s])]
        fullTrial_beam=[t for t, arr in enumerate(alignedData['beamAlign_isos_dFoF'][s]) if len(arr) == max(len(a) for a in alignedData['beamAlign_isos_dFoF'][s])]
        
        alignedData['cueAlign_isos_dFoF'][s] = ([alignedData['cueAlign_isos_dFoF'][s][t] for t in fullTrial_cue])
        alignedData['cueAlign_ACh_dFoF'][s] = ([alignedData['cueAlign_ACh_dFoF'][s][t] for t in fullTrial_cue])
        alignedData['cueAlign_DA_dFoF'][s] = ([alignedData['cueAlign_DA_dFoF'][s][t] for t in fullTrial_cue])
        alignedData['beamAlign_isos_dFoF'][s] = ([alignedData['beamAlign_isos_dFoF'][s][t] for t in fullTrial_beam])
        alignedData['beamAlign_ACh_dFoF'][s] = ([alignedData['beamAlign_ACh_dFoF'][s][t] for t in fullTrial_beam])
        alignedData['beamAlign_DA_dFoF'][s] = ([alignedData['beamAlign_DA_dFoF'][s][t] for t in fullTrial_beam])

    #%% Part 1.1.1 - remove trials with extreme movement artifacts

    threshold=8; alignmentOutput['index_trialWithArtifact']={}; alignmentOutput['index_trialWithoutArtifact']={}; alignmentOutput['index_trialWithArtifactResponses']={}; alignmentOutput['index_trialWithoutArtifactResponses']={}
    for s in range(len(sites)):
        alignmentOutput['index_trialWithArtifact'][s] = lf.findOutlierIndices(alignedData['cueAlign_ACh_dFoF'][s], threshold)
        
        if len(alignmentOutput['index_trialWithArtifact'][s]) == 1:
            fig,axs=plt.subplots(len(alignmentOutput['index_trialWithArtifact'][s]),1)
            plt.title('Identified Outliers')
            for trial in range(len(alignmentOutput['index_trialWithArtifact'][s])):
                plot1=axs.plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_isos_dFoF'][s][alignmentOutput['index_trialWithArtifact'][s][trial]], label='Trial {}'.format(alignmentOutput['index_trialWithArtifact'][s][trial]), color = 'k')
                #plot2=axs.plot(np.where(lf.identifyOutliers(alignedData['cueAlign_isos_dFoF'][alignmentOutput['index_trialWithArtifact'][s][trial]], threshold), trial, np.nan), 'ro')  # Plot outliers as red dots
                axs.legend()
                axs.set_xlabel('time (seconds)')
                axs.set_ylabel('isos dFoF')
            plt.show()    
        elif len(alignmentOutput['index_trialWithArtifact'][s]) > 0:
            fig,axs=plt.subplots(len(alignmentOutput['index_trialWithArtifact'][s]),1)
            plt.title('Identified Outliers')
            for trial in range(len(alignmentOutput['index_trialWithArtifact'][s])):
                plot1=axs[trial].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_isos_dFoF'][s][alignmentOutput['index_trialWithArtifact'][s][trial]], label='Trial {}'.format(alignmentOutput['index_trialWithArtifact'][s][trial]), color = 'k')
                #plot2=axs.plot(np.where(lf.identifyOutliers(alignedData['cueAlign_isos_dFoF'][alignmentOutput['index_trialWithArtifact'][s][trial]], threshold), trial, np.nan), 'ro')  # Plot outliers as red dots
                axs[trial].legend()
                axs[trial].set_xlabel('time (seconds)')
                axs[trial].set_ylabel('isos dFoF')
            plt.show()
          
        alignmentOutput['index_trialWithoutArtifact'][s] = [t for t in fullTrial_cue if t not in alignmentOutput['index_trialWithArtifact'][s]]
        artifactResponseIndices = []
        if (np.isin(np.where(alignmentOutput['index_trialWithResponse'][s]), alignmentOutput['index_trialWithArtifact'][s]).flatten()).any():
            artifactResponseIndices = np.argwhere(np.isin(np.where(alignmentOutput['index_trialWithResponse'][s]), alignmentOutput['index_trialWithArtifact'][s]).flatten() == True).flatten().tolist()
        
        alignmentOutput['index_trialWithArtifactResponses'][s] = artifactResponseIndices
        alignmentOutput['index_trialWithoutArtifactResponses'][s] = [t for t in alignmentOutput['index_trialWithResponse'][s] if t not in alignmentOutput['index_trialWithArtifactResponses'][s]]

        # if (np.isin(np.where(alignmentOutput['index_trialWithResponse']), alignmentOutput['index_trialWithArtifact'][s])).any():
        #     alignmentOutput['index_trialWithArtifactResponses'][s] = [np.where(np.isin(np.where(alignmentOutput['index_trialWithResponse']), alignmentOutput['index_trialWithArtifact'][s])==True)[1][0]]
        # else: # s for sites | t for trials
        #     alignmentOutput['index_trialWithArtifactResponses'][s] = []
            
        alignmentOutput['index_trialWithResponse'][s] = [alignmentOutput['index_trialWithResponse'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignmentOutput['index_trialWithStartCue'][s] = [alignmentOutput['index_trialWithStartCue'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignmentOutput['index_trialWithTone'][s] = [alignmentOutput['index_trialWithTone'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignmentOutput['index_trialWithTooEarlyResponse'][s] = [alignmentOutput['index_trialWithTooEarlyResponse'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignedData['cueAlign_isos_dFoF'][s] = [alignedData['cueAlign_isos_dFoF'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignedData['cueAlign_ACh_dFoF'][s] = [alignedData['cueAlign_ACh_dFoF'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignedData['cueAlign_DA_dFoF'][s] = [alignedData['cueAlign_DA_dFoF'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
        alignedData['beamAlign_isos_dFoF'][s] = [t for t, include in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]
        alignedData['beamAlign_ACh_dFoF'][s] = [t for t, include in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]
        alignedData['beamAlign_DA_dFoF'][s] = [t for t, include in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]

        if bxDataSheet['phaseIndex'][dayIndex] > 1:
            alignmentOutput['trialOutcome'][s] = [alignmentOutput['trialOutcome'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
            alignmentOutput['trialType'][s] = [alignmentOutput['trialType'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
            alignmentOutput['trialType_calc'][s] = [alignmentOutput['trialType_calc'][s][t] for t in alignmentOutput['index_trialWithoutArtifact'][s]]
            alignmentOutput['beamAlign_trialOutcome'][s] = [t for t, include in zip(alignmentOutput['beamAlign_trialOutcome'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]
            alignmentOutput['beamAlign_trialType'][s] = [t for t, include in zip(alignmentOutput['beamAlign_trialType'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]
            alignmentOutput['beamAlign_trialType_calc'][s] = [t for t, include in zip(alignmentOutput['beamAlign_trialType_calc'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include]
        
        if len(alignmentOutput['index_trialWithArtifactResponses'][s]) == 1:
            fig,axs=plt.subplots(len(alignmentOutput['index_trialWithArtifactResponses'][s]),1)
            plt.title('Identified Outliers')
            for trial in range(len(alignmentOutput['index_trialWithArtifactResponses'][s])): # s for sites
                plot1=axs.plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_isos_dFoF'][s][np.argwhere(np.isin(np.where(alignmentOutput['index_trialWithResponse'][s]),alignmentOutput['index_trialWithArtifactResponses'][s][trial]).flatten() == True).flatten().item()], label='Trial {}'.format(np.where(alignmentOutput['index_trialWithResponse'][s])[0][alignmentOutput['index_trialWithArtifactResponses'][s]][trial]), color = 'k')
                #plot2=axs[trial].plot(np.where(lf.identifyOutliers(alignedData['cueAlign_isos_dFoF'][alignmentOutput['index_trialWithArtifact'][s][trial]], threshold), trial, np.nan), 'ro')  # Plot outliers as red dots
                axs.legend()
                axs.set_xlabel('time (seconds)')
                axs.set_ylabel('isos dFoF')
            plt.show() 
        elif len(alignmentOutput['index_trialWithArtifactResponses'][s]) > 1:
            fig,axs=plt.subplots(len(alignmentOutput['index_trialWithArtifactResponses'][s]),1)
            plt.title('Identified Outliers')
            for trial in range(len(alignmentOutput['index_trialWithArtifactResponses'][s])): # s for sites
                plot1=axs[trial].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_isos_dFoF'][s][np.argwhere(np.isin(np.where(alignmentOutput['index_trialWithResponse'][s]),alignmentOutput['index_trialWithArtifactResponses'][s][trial]).flatten() == True).flatten().item()], label='Trial {}'.format(np.where(alignmentOutput['index_trialWithResponse'][s])[0][alignmentOutput['index_trialWithArtifactResponses'][s]][trial]), color = 'k')
                #plot2=axs[trial].plot(np.where(lf.identifyOutliers(alignedData['cueAlign_isos_dFoF'][alignmentOutput['index_trialWithArtifact'][s][trial]], threshold), trial, np.nan), 'ro')  # Plot outliers as red dots
                axs[trial].legend()
                axs[trial].set_xlabel('time (seconds)')
                axs[trial].set_ylabel('isos dFoF')
            plt.show()         
    
        alignedData['cueAlign_isos_stack'][s] = np.vstack(alignedData['cueAlign_isos_dFoF'][s])
        alignedData['cueAlign_ACh_stack'][s] = np.vstack(alignedData['cueAlign_ACh_dFoF'][s])
        alignedData['cueAlign_DA_stack'][s] = np.vstack(alignedData['cueAlign_DA_dFoF'][s])
        alignedData['beamAlign_isos_stack'][s] = np.vstack([t for t, include in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include])
        alignedData['beamAlign_ACh_stack'][s] = np.vstack([t for t, include in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include])
        alignedData['beamAlign_DA_stack'][s] = np.vstack([t for t, include in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['index_trialWithoutArtifactResponses'][s]) if include])

    #%% Part 1.1.2 - format data to cue-alignment    
    alignedData['cueAlign_isos_mean'] = {}; alignedData['cueAlign_isos_std'] = {}; alignedData['cueAlign_isos_mean_tooEarly'] = {}; alignedData['cueAlign_isos_std_tooEarly'] = {}; alignedData['cueAlign_isos_mean_goToneOnly'] = {}; alignedData['cueAlign_isos_std_goToneOnly'] = {}; alignedData['cueAlign_isos_mean_nogoToneOnly'] = {}; alignedData['cueAlign_isos_std_nogoToneOnly'] = {};
    alignedData['cueAlign_ACh_mean'] = {}; alignedData['cueAlign_ACh_std'] = {}; alignedData['cueAlign_ACh_mean_tooEarly'] = {}; alignedData['cueAlign_ACh_std_tooEarly'] = {}; alignedData['cueAlign_ACh_mean_goToneOnly'] = {}; alignedData['cueAlign_ACh_std_goToneOnly'] = {}; alignedData['cueAlign_ACh_mean_nogoToneOnly'] = {}; alignedData['cueAlign_ACh_std_nogoToneOnly'] = {};
    alignedData['cueAlign_DA_mean'] = {}; alignedData['cueAlign_DA_std'] = {}; alignedData['cueAlign_DA_mean_tooEarly'] = {}; alignedData['cueAlign_DA_std_tooEarly'] = {}; alignedData['cueAlign_DA_mean_goToneOnly'] = {}; alignedData['cueAlign_DA_std_goToneOnly'] = {}; alignedData['cueAlign_DA_mean_nogoToneOnly'] = {}; alignedData['cueAlign_DA_std_nogoToneOnly'] = {};
    alignedData['beamAlign_isos_mean'] = {}; alignedData['beamAlign_isos_std'] = {}; alignedData['beamAlign_isos_mean_goToneOnly'] = {}; alignedData['beamAlign_isos_std_goToneOnly'] = {}; alignedData['beamAlign_isos_mean_nogoToneOnly'] = {}; alignedData['beamAlign_isos_std_nogoToneOnly'] = {};
    alignedData['beamAlign_ACh_mean'] = {}; alignedData['beamAlign_ACh_std'] = {}; alignedData['beamAlign_ACh_mean_tooEarly'] = {}; alignedData['beamAlign_ACh_std_tooEarly'] = {}; alignedData['beamAlign_ACh_mean_goToneOnly'] = {}; alignedData['beamAlign_ACh_std_goToneOnly'] = {}; alignedData['beamAlign_ACh_mean_nogoToneOnly'] = {}; alignedData['beamAlign_ACh_std_nogoToneOnly'] = {};
    alignedData['beamAlign_DA_mean'] = {}; alignedData['beamAlign_DA_std'] = {}; alignedData['beamAlign_DA_mean_tooEarly'] = {}; alignedData['beamAlign_DA_std_tooEarly'] = {}; alignedData['beamAlign_DA_mean_goToneOnly'] = {}; alignedData['beamAlign_DA_std_goToneOnly'] = {}; alignedData['beamAlign_DA_mean_nogoToneOnly'] = {}; alignedData['beamAlign_DA_std_nogoToneOnly'] = {};
    for s in range(len(sites)):
        alignedData['cueAlign_isos_mean'][s] = np.mean(alignedData['cueAlign_isos_dFoF'][s], axis=0)
        alignedData['cueAlign_isos_std'][s] = np.std(alignedData['cueAlign_isos_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1: # go=1 nogo=0 early=2
            alignedData['cueAlign_isos_mean_tooEarly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0) 
            alignedData['cueAlign_isos_std_tooEarly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0)
            alignedData['cueAlign_isos_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['cueAlign_isos_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['cueAlign_isos_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['cueAlign_isos_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
        alignedData['cueAlign_ACh_mean'][s] = np.mean(alignedData['cueAlign_ACh_dFoF'][s], axis=0)
        alignedData['cueAlign_ACh_std'][s] = np.std(alignedData['cueAlign_ACh_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1:
            alignedData['cueAlign_ACh_mean_tooEarly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0)
            alignedData['cueAlign_ACh_std_tooEarly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0)
            alignedData['cueAlign_ACh_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['cueAlign_ACh_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['cueAlign_ACh_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['cueAlign_ACh_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
        alignedData['cueAlign_DA_mean'][s] = np.mean(alignedData['cueAlign_DA_dFoF'][s], axis=0)
        alignedData['cueAlign_DA_std'][s] = np.std(alignedData['cueAlign_DA_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1:
            alignedData['cueAlign_DA_mean_tooEarly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0)
            alignedData['cueAlign_DA_std_tooEarly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 2], axis=0)
            alignedData['cueAlign_DA_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['cueAlign_DA_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['cueAlign_DA_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['cueAlign_DA_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s],alignmentOutput['trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
        alignedData['beamAlign_isos_mean'][s] = np.mean(alignedData['beamAlign_isos_dFoF'][s], axis=0)
        alignedData['beamAlign_isos_std'][s] = np.std(alignedData['beamAlign_isos_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1:
            alignedData['beamAlign_isos_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['beamAlign_isos_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['beamAlign_isos_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['beamAlign_isos_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
        alignedData['beamAlign_ACh_mean'][s] = np.mean(alignedData['beamAlign_ACh_dFoF'][s], axis=0)
        alignedData['beamAlign_ACh_std'][s] = np.std(alignedData['beamAlign_ACh_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1:
            alignedData['beamAlign_ACh_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['beamAlign_ACh_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['beamAlign_ACh_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['beamAlign_ACh_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
        alignedData['beamAlign_DA_mean'][s] = np.mean(alignedData['beamAlign_DA_dFoF'][s], axis=0)
        alignedData['beamAlign_DA_std'][s] = np.std(alignedData['beamAlign_DA_dFoF'][s], axis=0)
        if keyDataSheet['exp_seshPhase'] > 1:
            alignedData['beamAlign_DA_mean_goToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
            alignedData['beamAlign_DA_std_goToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 1], axis=0)
        if keyDataSheet['exp_seshPhase'] > 2:
            alignedData['beamAlign_DA_mean_nogoToneOnly'][s] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
            alignedData['beamAlign_DA_std_nogoToneOnly'][s] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s],alignmentOutput['beamAlign_trialType_calc'][s]) if thisTrialType == 0], axis=0)
        
    #%% Part 1.1.3 - format data to BINNED cue-alignment
    
    nBins = 3
    trialBins = int(np.ceil(np.size(alignedData['cueAlign_isos_dFoF'][0],0)/nBins))
    trialBinsBeam = int(np.ceil(np.size(alignedData['beamAlign_isos_dFoF'][0],0)/nBins))
    
    alignmentOutputBinned = {'cueAlign_isos_mean' : {}, 'cueAlign_isos_std' : {}, 'cueAlign_isos_mean_tooEarly' : {}, 'cueAlign_isos_std_tooEarly' : {}, 'cueAlign_isos_mean_goToneOnly' : {}, 'cueAlign_isos_std_goToneOnly' : {}, 'cueAlign_isos_mean_nogoToneOnly' : {}, 'cueAlign_isos_std_nogoToneOnly' : {}, 'cueAlign_ACh_mean' : {}, 'cueAlign_ACh_std' : {}, 'cueAlign_ACh_mean_tooEarly' : {}, 'cueAlign_ACh_std_tooEarly' : {}, 'cueAlign_ACh_mean_goToneOnly' : {}, 'cueAlign_ACh_std_goToneOnly' : {}, 'cueAlign_ACh_mean_nogoToneOnly' : {}, 'cueAlign_ACh_std_nogoToneOnly' : {}, 'cueAlign_DA_mean' : {}, 'cueAlign_DA_std' : {}, 'cueAlign_DA_mean_tooEarly' : {}, 'cueAlign_DA_std_tooEarly' : {}, 'cueAlign_DA_mean_goToneOnly' : {}, 'cueAlign_DA_std_goToneOnly' : {}, 'cueAlign_DA_mean_nogoToneOnly' : {}, 'cueAlign_DA_std_nogoToneOnly' : {}, 'beamAlign_isos_mean' : {}, 'beamAlign_isos_std' : {}, 'beamAlign_isos_mean_tooEarly' : {}, 'beamAlign_isos_std_tooEarly' : {}, 'beamAlign_isos_mean_goToneOnly' : {}, 'beamAlign_isos_std_goToneOnly' : {}, 'beamAlign_isos_mean_nogoToneOnly' : {}, 'beamAlign_isos_std_nogoToneOnly' : {}, 'beamAlign_ACh_mean' : {}, 'beamAlign_ACh_std' : {}, 'beamAlign_ACh_mean_tooEarly' : {}, 'beamAlign_ACh_std_tooEarly' : {}, 'beamAlign_ACh_mean_goToneOnly' : {}, 'beamAlign_ACh_std_goToneOnly' : {}, 'beamAlign_ACh_mean_nogoToneOnly' : {}, 'beamAlign_ACh_std_nogoToneOnly' : {}, 'beamAlign_DA_mean' : {}, 'beamAlign_DA_std' : {}, 'beamAlign_DA_mean_tooEarly' : {}, 'beamAlign_DA_std_tooEarly' : {}, 'beamAlign_DA_mean_goToneOnly' : {}, 'beamAlign_DA_std_goToneOnly' : {}, 'beamAlign_DA_mean_nogoToneOnly' : {}, 'beamAlign_DA_std_nogoToneOnly' : {}}
    for s in range(len(sites)):
        tTrialsInBin = 0
        tBeamTrialsInBin = 0
        alignmentOutputBinned['cueAlign_isos_mean'][s] = {}; alignmentOutputBinned['cueAlign_isos_std'][s] = {}; alignmentOutputBinned['cueAlign_isos_mean_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_isos_std_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_isos_mean_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_isos_std_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_isos_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_isos_std_nogoToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_mean'][s] = {}; alignmentOutputBinned['cueAlign_ACh_std'][s] = {}; alignmentOutputBinned['cueAlign_ACh_mean_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_std_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_mean_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_std_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_ACh_std_nogoToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_DA_mean'][s] = {}; alignmentOutputBinned['cueAlign_DA_std'][s] = {}; alignmentOutputBinned['cueAlign_DA_mean_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_DA_std_tooEarly'][s] = {}; alignmentOutputBinned['cueAlign_DA_mean_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_DA_std_goToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_DA_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['cueAlign_DA_std_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_isos_mean'][s] = {}; alignmentOutputBinned['beamAlign_isos_std'][s] = {}; alignmentOutputBinned['beamAlign_isos_mean_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_isos_std_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_isos_mean_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_isos_std_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_isos_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_isos_std_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_mean'][s] = {}; alignmentOutputBinned['beamAlign_ACh_std'][s] = {}; alignmentOutputBinned['beamAlign_ACh_mean_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_std_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_mean_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_std_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_ACh_std_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_DA_mean'][s] = {}; alignmentOutputBinned['beamAlign_DA_std'][s] = {}; alignmentOutputBinned['beamAlign_DA_mean_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_DA_std_tooEarly'][s] = {}; alignmentOutputBinned['beamAlign_DA_mean_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_DA_std_goToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_DA_mean_nogoToneOnly'][s] = {}; alignmentOutputBinned['beamAlign_DA_std_nogoToneOnly'][s] = {};        
        for binIndex in range(nBins):
            alignmentOutputBinned['cueAlign_isos_mean'][s][binIndex] = np.mean(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            alignmentOutputBinned['cueAlign_isos_std'][s][binIndex] = np.std(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['cueAlign_isos_mean_tooEarly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_isos_std_tooEarly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_isos_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['cueAlign_isos_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['cueAlign_isos_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['cueAlign_isos_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_isos_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
            
            alignmentOutputBinned['cueAlign_ACh_mean'][s][binIndex] = np.mean(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            alignmentOutputBinned['cueAlign_ACh_std'][s][binIndex] = np.std(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['cueAlign_ACh_mean_tooEarly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_ACh_std_tooEarly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_ACh_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['cueAlign_ACh_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['cueAlign_ACh_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['cueAlign_ACh_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_ACh_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
            
            alignmentOutputBinned['cueAlign_DA_mean'][s][binIndex] = np.mean(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            alignmentOutputBinned['cueAlign_DA_std'][s][binIndex] = np.std(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['cueAlign_DA_mean_tooEarly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_DA_std_tooEarly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 2], axis=0)
                alignmentOutputBinned['cueAlign_DA_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['cueAlign_DA_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['cueAlign_DA_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['cueAlign_DA_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['cueAlign_DA_dFoF'][s][tTrialsInBin:tTrialsInBin+trialBins],alignmentOutput['trialType_calc'][s][tTrialsInBin:tTrialsInBin+trialBins]) if thisTrialType == 0], axis=0)
            
            alignmentOutputBinned['beamAlign_isos_mean'][s][binIndex] = np.mean(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            alignmentOutputBinned['beamAlign_isos_std'][s][binIndex] = np.std(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['beamAlign_isos_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['beamAlign_isos_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['beamAlign_isos_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['beamAlign_isos_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_isos_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
            
            alignmentOutputBinned['beamAlign_ACh_mean'][s][binIndex] = np.mean(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            alignmentOutputBinned['beamAlign_ACh_std'][s][binIndex] = np.std(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['beamAlign_ACh_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['beamAlign_ACh_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['beamAlign_ACh_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['beamAlign_ACh_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_ACh_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
            
            alignmentOutputBinned['beamAlign_DA_mean'][s][binIndex] = np.mean(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            alignmentOutputBinned['beamAlign_DA_std'][s][binIndex] = np.std(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam], axis=0)
            if keyDataSheet['exp_seshPhase'] > 1:
                alignmentOutputBinned['beamAlign_DA_mean_goToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
                alignmentOutputBinned['beamAlign_DA_std_goToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 1], axis=0)
            if keyDataSheet['exp_seshPhase'] > 2:
                alignmentOutputBinned['beamAlign_DA_mean_nogoToneOnly'][s][binIndex] = np.mean([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
                alignmentOutputBinned['beamAlign_DA_std_nogoToneOnly'][s][binIndex] = np.std([t for t, thisTrialType in zip(alignedData['beamAlign_DA_dFoF'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam],alignmentOutput['beamAlign_trialType_calc'][s][tBeamTrialsInBin:tBeamTrialsInBin+trialBinsBeam]) if thisTrialType == 0], axis=0)
        
            #print(list(range(tTrialsInBin,tTrialsInBin+trialBins)))
            tTrialsInBin = tTrialsInBin + trialBins
            #print(list(range(tBeamTrialsInBin,tBeamTrialsInBin+trialBinsBeam)))
            tBeamTrialsInBin = tBeamTrialsInBin + trialBinsBeam


    #%% Part 1.2 - plotting cue-alignment
    for s in range(len(sites)):
        print('    Plotting {} hemis cue align PSTH:'.format(sites[s]))
        
        ### !! GLOBAL CUE !!
        fig,axs=plt.subplots(2,1)
        vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
        if keyDataSheet['exp_seshPhase'] != 1:
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
        plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean'][s], 'r', label='cue-align DA')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean'][s] - alignedData['cueAlign_DA_std'][s], alignedData['cueAlign_DA_mean'][s] + alignedData['cueAlign_DA_std'][s], color='r', alpha=0.3)
        
        axs[0].set_xlabel('time from start (seconds)')
        axs[0].set_ylabel('dF/F (%)')
        axs[0].set_title('DA dF/F - cue align - {}'.format(sites[s]))# - {} trials'.format())
        if keyDataSheet['exp_seshPhase'] == 1:
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
        else:
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                        Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
        axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[0].set_xlim(-2,6)
        axs[0].set_ylim(-2,3)
        
        vline1=axs[1].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
        if keyDataSheet['exp_seshPhase'] != 1:
            vline2=axs[1].axvline(x=2.5, color='k', linestyle='--', label='tone')
        plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean'][s], 'g', label='cue-align ACh')
        axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean'][s] - alignedData['cueAlign_ACh_std'][s], alignedData['cueAlign_ACh_mean'][s] + alignedData['cueAlign_ACh_std'][s], color='g', alpha=0.3)
        
        axs[1].set_xlabel('time from start (seconds)')
        axs[1].set_ylabel('dF/F (%)')
        axs[1].set_title('ACh dF/F - cue align - {}'.format(sites[s]))
        if keyDataSheet['exp_seshPhase'] == 1:
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
        else:
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                        Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
        axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[1].set_xlim(-2,6)
        axs[1].set_ylim(-2,3)
        
        plt.savefig('{directory}cueAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        ### !! GLOBAL BEAM !!
        fig,axs=plt.subplots(2,1)
        vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
        plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean'][s], 'r', label='beam-align DA')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean'][s] - alignedData['beamAlign_DA_std'][s], alignedData['beamAlign_DA_mean'][s] + alignedData['beamAlign_DA_std'][s], color='r', alpha=0.3)
        
        axs[0].set_xlabel('time from start (seconds)')
        axs[0].set_ylabel('dF/F (%)')
        axs[0].set_title('DA dF/F - beam break - {}'.format(sites[s]))# - {} trials'.format())
        lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                        Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
        axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[0].set_xlim(-2,3)
        axs[0].set_ylim(-2,3)
        
        vline1=axs[1].axvline(x=0, color='b', linestyle='--', label='beam break')
        plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean'][s], 'g', label='beam-align ACh')
        axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean'][s] - alignedData['beamAlign_ACh_std'][s], alignedData['beamAlign_ACh_mean'][s] + alignedData['beamAlign_ACh_std'][s], color='g', alpha=0.3)
        
        axs[1].set_xlabel('time from start (seconds)')
        axs[1].set_ylabel('dF/F (%)')
        axs[1].set_title('ACh dF/F - beam break - {}'.format(sites[s]))
        lines = [Line2D([0], [0], color='g', lw=2, label='beam-align ACh'),
                        Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
        axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[1].set_xlim(-2,3)
        axs[1].set_ylim(-2,3)
        
        plt.savefig('{directory}beamAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        if keyDataSheet['exp_seshPhase'] > 1:
            ### !! TOO EARLY !!
            fig,axs=plt.subplots(2,1)
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
    #        vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_tooEarly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_tooEarly'][s] - alignedData['cueAlign_DA_std_tooEarly'][s], alignedData['cueAlign_DA_mean_tooEarly'][s] + alignedData['cueAlign_DA_std_tooEarly'][s], color='r', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('DA dF/F - too Early - {}'.format(sites[s]))# - {} trials'.format())
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            vline1=axs[1].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
    #        vline2=axs[1].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_tooEarly'][s], 'g', label='cue-align ACh')
            axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_tooEarly'][s] - alignedData['cueAlign_ACh_std_tooEarly'][s], alignedData['cueAlign_ACh_mean_tooEarly'][s] + alignedData['cueAlign_ACh_std_tooEarly'][s], color='g', alpha=0.3)
        
            axs[1].set_xlabel('time from start (seconds)')
            axs[1].set_ylabel('dF/F (%)')
            axs[1].set_title('ACh dF/F - too Early - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
            axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[1].set_xlim(-2,6)
            axs[1].set_ylim(-2,3)
        
            plt.savefig('{directory}cueAlign_tooEarly_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
        ### !! GO TONE !!
        if keyDataSheet['exp_seshPhase'] > 1:
            fig,axs=plt.subplots(2,1)
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_goToneOnly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_goToneOnly'][s] - alignedData['cueAlign_DA_std_goToneOnly'][s], alignedData['cueAlign_DA_mean_goToneOnly'][s] + alignedData['cueAlign_DA_std_goToneOnly'][s], color='r', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('DA dF/F - go Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            vline1=axs[1].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[1].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_goToneOnly'][s], 'g', label='cue-align ACh')
            axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_goToneOnly'][s] - alignedData['cueAlign_ACh_std_goToneOnly'][s], alignedData['cueAlign_ACh_mean_goToneOnly'][s] + alignedData['cueAlign_ACh_std_goToneOnly'][s], color='g', alpha=0.3)
        
            axs[1].set_xlabel('time from start (seconds)')
            axs[1].set_ylabel('dF/F (%)')
            axs[1].set_title('ACh dF/F - go Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[1].set_xlim(-2,6)
            axs[1].set_ylim(-2,3)
        
            plt.savefig('{directory}cueAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM 
            fig,axs=plt.subplots(2,1)
            vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_goToneOnly'][s], 'r', label='beam-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_goToneOnly'][s] - alignedData['beamAlign_DA_std_goToneOnly'][s], alignedData['beamAlign_DA_mean_goToneOnly'][s] + alignedData['beamAlign_DA_std_goToneOnly'][s], color='r', alpha=0.3)
        
            axs[0].set_xlabel('time from beam break (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('DA dF/F - Beam Break After go Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,3)
            axs[0].set_ylim(-2,3)
        
            vline1=axs[1].axvline(x=0, color='b', linestyle='--', label='beam break')
            plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_goToneOnly'][s], 'g', label='beam-align ACh')
            axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_goToneOnly'][s] - alignedData['beamAlign_ACh_std_goToneOnly'][s], alignedData['beamAlign_ACh_mean_goToneOnly'][s] + alignedData['beamAlign_ACh_std_goToneOnly'][s], color='g', alpha=0.3)
        
            axs[1].set_xlabel('time from beam break (seconds)')
            axs[1].set_ylabel('dF/F (%)')
            axs[1].set_title('ACh dF/F - Beam Break After go Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[1].set_xlim(-2,3)
            axs[1].set_ylim(-2,3)
        
            plt.savefig('{directory}beamAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! NOGO TONE !!
        if keyDataSheet['exp_seshPhase'] > 2:
            fig,axs=plt.subplots(2,1)
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_nogoToneOnly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_nogoToneOnly'][s] - alignedData['cueAlign_DA_std_nogoToneOnly'][s], alignedData['cueAlign_DA_mean_nogoToneOnly'][s] + alignedData['cueAlign_DA_std_nogoToneOnly'][s], color='r', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('DA dF/F - nogo Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            vline1=axs[1].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[1].axvline(x=2.5, color='k', linestyle='--', label='tone')
            plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s], 'g', label='cue-align ACh')
            axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s] - alignedData['cueAlign_ACh_std_nogoToneOnly'][s], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s] + alignedData['cueAlign_ACh_std_nogoToneOnly'][s], color='g', alpha=0.3)
        
            axs[1].set_xlabel('time from start (seconds)')
            axs[1].set_ylabel('dF/F (%)')
            axs[1].set_title('ACh dF/F - nogo Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[1].set_xlim(-2,6)
            axs[1].set_ylim(-2,3)
        
            plt.savefig('{directory}cueAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM 
            fig,axs=plt.subplots(2,1)
            vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_nogoToneOnly'][s], 'r', label='beam-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_nogoToneOnly'][s] - alignedData['beamAlign_DA_std_nogoToneOnly'][s], alignedData['beamAlign_DA_mean_nogoToneOnly'][s] + alignedData['beamAlign_DA_std_nogoToneOnly'][s], color='r', alpha=0.3)
        
            axs[0].set_xlabel('time from beam break (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('DA dF/F - Beam Break After nogo Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,3)
            axs[0].set_ylim(-2,3)
        
            vline1=axs[1].axvline(x=0, color='b', linestyle='--', label='beam break')
            plot2=axs[1].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s], 'g', label='beam-align ACh')
            axs[1].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s] - alignedData['beamAlign_ACh_std_nogoToneOnly'][s], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s] + alignedData['beamAlign_ACh_std_nogoToneOnly'][s], color='g', alpha=0.3)
        
            axs[1].set_xlabel('time from beam break (seconds)')
            axs[1].set_ylabel('dF/F (%)')
            axs[1].set_title('ACh dF/F - Beam Break After nogo Tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[1].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[1].set_xlim(-2,3)
            axs[1].set_ylim(-2,3)
        
            plt.savefig('{directory}beamAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            

    #%% Part 1.3 - plotting STACKED cue-alignment
    for s in range(len(sites)):
        print('    Plotting {} hemis stacked align PSTH:'.format(sites[s]))
        
        ### !! GLOBAL CUE !!
        fig,axs=plt.subplots(1,1); axs = [axs]
        vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
        if keyDataSheet['exp_seshPhase'] != 1:
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
        
        plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean'][s], 'r', label='cue-align DA')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean'][s] - alignedData['cueAlign_DA_std'][s], alignedData['cueAlign_DA_mean'][s] + alignedData['cueAlign_DA_std'][s], color='r', alpha=0.3)
        plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean'][s], 'g', label='cue-align ACh')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean'][s] - alignedData['cueAlign_ACh_std'][s], alignedData['cueAlign_ACh_mean'][s] + alignedData['cueAlign_ACh_std'][s], color='g', alpha=0.3)
        #plot3=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_isos_mean'][s], 'b', label='cue-align isos')
        #axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_isos_mean'][s] - alignedData['cueAlign_isos_std'][s], alignedData['cueAlign_isos_mean'][s] + alignedData['cueAlign_isos_std'][s], color='b', alpha=0.3)
        
        axs[0].set_xlabel('time from start (seconds)')
        axs[0].set_ylabel('dF/F (%)')
        axs[0].set_title('cue align - {}'.format(sites[s]))# - {} trials'.format())
        if keyDataSheet['exp_seshPhase'] == 1:
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
        #             Line2D([0], [0], color='b', lw=2, label='cue-align isos'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
        else:
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
        #             Line2D([0], [0], color='b', lw=2, label='cue-align isos'),
                        Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                        Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
        axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[0].set_xlim(-2,6)
        axs[0].set_ylim(-2,3)
        
        plt.savefig('{directory}stackedCueAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        ### !! GLOBAL BEAM !!
        fig,axs=plt.subplots(1,1); axs = [axs]
        vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
        
        plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean'][s], 'r', label='beam-align DA')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean'][s] - alignedData['beamAlign_DA_std'][s], alignedData['beamAlign_DA_mean'][s] + alignedData['beamAlign_DA_std'][s], color='r', alpha=0.3)
        plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean'][s], 'g', label='beam-align ACh')
        axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean'][s] - alignedData['beamAlign_ACh_std'][s], alignedData['beamAlign_ACh_mean'][s] + alignedData['beamAlign_ACh_std'][s], color='g', alpha=0.3)
        
        axs[0].set_xlabel('time from start (seconds)')
        axs[0].set_ylabel('dF/F (%)')
        axs[0].set_title('beam break align - {}'.format(sites[s]))# - {} trials'.format())
        lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                 Line2D([0], [0], color='g', lw=2, label='beam-align ACh'),
                        Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
        axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
        axs[0].set_xlim(-2,3)
        axs[0].set_ylim(-2,3)
        
        plt.savefig('{directory}stackedBeamAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        if keyDataSheet['exp_seshPhase'] > 1:
            ### !! TOO EARLY !!
            fig,axs=plt.subplots(1,1); axs = [axs]
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
    #        vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_tooEarly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_tooEarly'][s] - alignedData['cueAlign_DA_std_tooEarly'][s], alignedData['cueAlign_DA_mean_tooEarly'][s] + alignedData['cueAlign_DA_std_tooEarly'][s], color='r', alpha=0.3)
            plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_tooEarly'][s], 'g', label='cue-align ACh')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_tooEarly'][s] - alignedData['cueAlign_ACh_std_tooEarly'][s], alignedData['cueAlign_ACh_mean_tooEarly'][s] + alignedData['cueAlign_ACh_std_tooEarly'][s], color='g', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('too Early - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            plt.savefig('{directory}stackedCueAlign_tooEarly_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
           
        ### IF TRAINING IS BEYOND REACHING PHASE    
        if keyDataSheet['exp_seshPhase'] > 1:
            ### !! GO TONE !!
            fig,axs=plt.subplots(1,1); axs = [axs]
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_goToneOnly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_goToneOnly'][s] - alignedData['cueAlign_DA_std_goToneOnly'][s], alignedData['cueAlign_DA_mean_goToneOnly'][s] + alignedData['cueAlign_DA_std_goToneOnly'][s], color='r', alpha=0.3)
            plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_goToneOnly'][s], 'g', label='cue-align ACh')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_goToneOnly'][s] - alignedData['cueAlign_ACh_std_goToneOnly'][s], alignedData['cueAlign_ACh_mean_goToneOnly'][s] + alignedData['cueAlign_ACh_std_goToneOnly'][s], color='g', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('go tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            plt.savefig('{directory}stackedCueAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM !!
            fig,axs=plt.subplots(1,1); axs = [axs]
            vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
            
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_goToneOnly'][s], 'r', label='beam-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_goToneOnly'][s] - alignedData['beamAlign_DA_std_goToneOnly'][s], alignedData['beamAlign_DA_mean_goToneOnly'][s] + alignedData['beamAlign_DA_std_goToneOnly'][s], color='r', alpha=0.3)
            plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_goToneOnly'][s], 'g', label='beam-align ACh')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_goToneOnly'][s] - alignedData['beamAlign_ACh_std_goToneOnly'][s], alignedData['beamAlign_ACh_mean_goToneOnly'][s] + alignedData['beamAlign_ACh_std_goToneOnly'][s], color='g', alpha=0.3)
        
            axs[0].set_xlabel('time from beam break (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('beam break after go tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='beam-align ACh'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,3)
            axs[0].set_ylim(-2,3)
        
            plt.savefig('{directory}stackedBeamAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            plt.close()
            
            ### !! NOGO TONE !!
        if keyDataSheet['exp_seshPhase'] > 2:
            fig,axs=plt.subplots(1,1); axs = [axs]
            vline1=axs[0].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            vline2=axs[0].axvline(x=2.5, color='k', linestyle='--', label='tone')
            
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_nogoToneOnly'][s], 'r', label='cue-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_DA_mean_nogoToneOnly'][s] - alignedData['cueAlign_DA_std_nogoToneOnly'][s], alignedData['cueAlign_DA_mean_nogoToneOnly'][s] + alignedData['cueAlign_DA_std_nogoToneOnly'][s], color='r', alpha=0.3)
            plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s], 'g', label='cue-align ACh')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s] - alignedData['cueAlign_ACh_std_nogoToneOnly'][s], alignedData['cueAlign_ACh_mean_nogoToneOnly'][s] + alignedData['cueAlign_ACh_std_nogoToneOnly'][s], color='g', alpha=0.3)
        
            axs[0].set_xlabel('time from start (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('nogo tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='cue-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([0], [0], color='k', linestyle='--', lw=2, label='tone')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,6)
            axs[0].set_ylim(-2,3)
        
            plt.savefig('{directory}stackedCueAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM !!
            fig,axs=plt.subplots(1,1); axs = [axs]
            vline1=axs[0].axvline(x=0, color='b', linestyle='--', label='beam break')
            
            plot1=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_nogoToneOnly'][s], 'r', label='beam-align DA')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_DA_mean_nogoToneOnly'][s] - alignedData['beamAlign_DA_std_nogoToneOnly'][s], alignedData['beamAlign_DA_mean_nogoToneOnly'][s] + alignedData['beamAlign_DA_std_nogoToneOnly'][s], color='r', alpha=0.3)
            plot2=axs[0].plot(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s], 'g', label='beam-align ACh')
            axs[0].fill_between(alignmentOutput['xAxisPoints'], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s] - alignedData['beamAlign_ACh_std_nogoToneOnly'][s], alignedData['beamAlign_ACh_mean_nogoToneOnly'][s] + alignedData['beamAlign_ACh_std_nogoToneOnly'][s], color='g', alpha=0.3)
        
            axs[0].set_xlabel('time from beam break (seconds)')
            axs[0].set_ylabel('dF/F (%)')
            axs[0].set_title('beam break after nogo tone - {}'.format(sites[s]))
            lines = [Line2D([0], [0], color='r', lw=2, label='beam-align DA'),
                     Line2D([0], [0], color='g', lw=2, label='cue-align ACh'),
                            Line2D([0], [0], color='b', linestyle='--', lw=2, label='beam break')]
            axs[0].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[0].set_xlim(-2,3)
            axs[0].set_ylim(-2,3)
        
            plt.savefig('{directory}stackedBeamAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
 
    #%% Part 1.4 - plotting BINNED cue-alignment
    for s in range(len(sites)):
        print('    Plotting {} hemis binned align PSTH:'.format(sites[s]))
   
        ### !! GLOBAL CUE !!
        fig,axs=plt.subplots(nBins,1)
        axs[0].set_title('cue align - {}'.format(sites[s]))# - {} trials'.format())
        for binIndex in range(nBins):
            vline1=axs[binIndex].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
            if keyDataSheet['exp_seshPhase'] != 1:
                vline2=axs[binIndex].axvline(x=2.5, color='k', linestyle='--', label='tone')
            
            plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean'][s][binIndex], 'r', label='cue-align DA')
            axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean'][s][binIndex] - alignmentOutputBinned['cueAlign_DA_std'][s][binIndex], alignmentOutputBinned['cueAlign_DA_mean'][s][binIndex] + alignmentOutputBinned['cueAlign_DA_std'][s][binIndex], color='r', alpha=0.3)
            plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean'][s][binIndex], 'g', label='cue-align ACh')
            axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean'][s][binIndex] - alignmentOutputBinned['cueAlign_ACh_std'][s][binIndex], alignmentOutputBinned['cueAlign_ACh_mean'][s][binIndex] + alignmentOutputBinned['cueAlign_ACh_std'][s][binIndex], color='g', alpha=0.3)
            
            axs[binIndex].set_xlabel('time from start (seconds)')
            axs[binIndex].set_ylabel('dF/F (%)')
            if keyDataSheet['exp_seshPhase'] == 1:
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='cue-align ACh'),
                            Line2D([binIndex], [binIndex], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
            else:
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='cue-align ACh'),
                            Line2D([binIndex], [binIndex], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                            Line2D([binIndex], [binIndex], color='k', linestyle='--', lw=2, label='tone')]
            axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[binIndex].set_xlim(-2,6)
            axs[binIndex].set_ylim(-2,3)
            
        plt.savefig('{directory}binnedCueAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        ### !! GLOBAL BEAM !!
        fig,axs=plt.subplots(nBins,1)
        axs[0].set_title('beam align - {}'.format(sites[s]))# - {} trials'.format())
        for binIndex in range(nBins):
            vline1=axs[binIndex].axvline(x=0, color='b', linestyle='--', label='beam break')
        
            plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean'][s][binIndex], 'r', label='beam-align DA')
            axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean'][s][binIndex] - alignmentOutputBinned['beamAlign_DA_std'][s][binIndex], alignmentOutputBinned['beamAlign_DA_mean'][s][binIndex] + alignmentOutputBinned['beamAlign_DA_std'][s][binIndex], color='r', alpha=0.3)
            plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean'][s][binIndex], 'g', label='beam-align ACh')
            axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean'][s][binIndex] - alignmentOutputBinned['beamAlign_ACh_std'][s][binIndex], alignmentOutputBinned['beamAlign_ACh_mean'][s][binIndex] + alignmentOutputBinned['beamAlign_ACh_std'][s][binIndex], color='g', alpha=0.3)
            
            axs[binIndex].set_xlabel('time from start (seconds)')
            axs[binIndex].set_ylabel('dF/F (%)')
            lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='beam-align DA'),
                     Line2D([binIndex], [binIndex], color='g', lw=2, label='beam-align ACh'),
                            Line2D([binIndex], [binIndex], color='b', linestyle='--', lw=2, label='beam break')]
            axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
            axs[binIndex].set_xlim(-2,3)
            axs[binIndex].set_ylim(-2,3)
            
        plt.savefig('{directory}binnedBeamAlign_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        #plt.close()
        
        if keyDataSheet['exp_seshPhase'] > 1:
            ### !! TOO EARLY !!
            fig,axs=plt.subplots(nBins,1)
            axs[0].set_title('too Early - {}'.format(sites[s]))
            for binIndex in range(nBins):
                vline1=axs[binIndex].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
                #vline2=axs[binIndex].axvline(x=2.5, color='k', linestyle='--', label='tone')
            
                plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_tooEarly'][s][binIndex], 'r', label='cue-align DA')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_tooEarly'][s][binIndex] - alignmentOutputBinned['cueAlign_DA_std_tooEarly'][s][binIndex], alignmentOutputBinned['cueAlign_DA_mean_tooEarly'][s][binIndex] + alignmentOutputBinned['cueAlign_DA_std_tooEarly'][s][binIndex], color='r', alpha=0.3)
                plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_tooEarly'][s][binIndex], 'g', label='cue-align ACh')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_tooEarly'][s][binIndex] - alignmentOutputBinned['cueAlign_ACh_std_tooEarly'][s][binIndex], alignmentOutputBinned['cueAlign_ACh_mean_tooEarly'][s][binIndex] + alignmentOutputBinned['cueAlign_ACh_std_tooEarly'][s][binIndex], color='g', alpha=0.3)
            
                axs[binIndex].set_xlabel('time from start (seconds)')
                axs[binIndex].set_ylabel('dF/F (%)')
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='cue-align ACh'),
                                Line2D([binIndex], [binIndex], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue')]
                axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
                axs[binIndex].set_xlim(-2,6)
                axs[binIndex].set_ylim(-2,3)
            
            plt.savefig('{directory}binnedCueAlign_tooEarly_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
           
        ### IF TRAINING IS BEYOND REACHING PHASE    
        if keyDataSheet['exp_seshPhase'] > 1:
            ### !! GO TONE BINNED !!
            fig,axs=plt.subplots(nBins,1)
            axs[0].set_title('go tone - {}'.format(sites[s]))
            for binIndex in range(nBins):
                vline1=axs[binIndex].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
                vline2=axs[binIndex].axvline(x=2.5, color='k', linestyle='--', label='tone')
                
                plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_goToneOnly'][s][binIndex], 'r', label='cue-align DA')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_goToneOnly'][s][binIndex] - alignmentOutputBinned['cueAlign_DA_std_goToneOnly'][s][binIndex], alignmentOutputBinned['cueAlign_DA_mean_goToneOnly'][s][binIndex] + alignmentOutputBinned['cueAlign_DA_std_goToneOnly'][s][binIndex], color='r', alpha=0.3)
                plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_goToneOnly'][s][binIndex], 'g', label='cue-align ACh')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_goToneOnly'][s][binIndex] - alignmentOutputBinned['cueAlign_ACh_std_goToneOnly'][s][binIndex], alignmentOutputBinned['cueAlign_ACh_mean_goToneOnly'][s][binIndex] + alignmentOutputBinned['cueAlign_ACh_std_goToneOnly'][s][binIndex], color='g', alpha=0.3)
            
                axs[binIndex].set_xlabel('time from start (seconds)')
                axs[binIndex].set_ylabel('dF/F (%)')
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='cue-align ACh'),
                                Line2D([binIndex], [binIndex], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                                Line2D([binIndex], [binIndex], color='k', linestyle='--', lw=2, label='tone')]
                axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
                axs[binIndex].set_xlim(-2,6)
                axs[binIndex].set_ylim(-2,3)
                        
            plt.savefig('{directory}binnedCueAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM 
            fig,axs=plt.subplots(nBins,1)
            axs[0].set_title('beam break after go tone - {}'.format(sites[s]))
            for binIndex in range(nBins):
                vline1=axs[binIndex].axvline(x=0, color='b', linestyle='--', label='beam break')
            
                plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean_goToneOnly'][s][binIndex], 'r', label='beam-align DA')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean_goToneOnly'][s][binIndex] - alignmentOutputBinned['beamAlign_DA_std_goToneOnly'][s][binIndex], alignmentOutputBinned['beamAlign_DA_mean_goToneOnly'][s][binIndex] + alignmentOutputBinned['beamAlign_DA_std_goToneOnly'][s][binIndex], color='r', alpha=0.3)
                plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean_goToneOnly'][s][binIndex], 'g', label='beam-align ACh')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean_goToneOnly'][s][binIndex] - alignmentOutputBinned['beamAlign_ACh_std_goToneOnly'][s][binIndex], alignmentOutputBinned['beamAlign_ACh_mean_goToneOnly'][s][binIndex] + alignmentOutputBinned['beamAlign_ACh_std_goToneOnly'][s][binIndex], color='g', alpha=0.3)
            
                axs[binIndex].set_xlabel('time from beam break (seconds)')
                axs[binIndex].set_ylabel('dF/F (%)')
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='beam-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='beam-align ACh'),
                                Line2D([binIndex], [binIndex], color='b', linestyle='--', lw=2, label='beam break')]
                axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
                axs[binIndex].set_xlim(-2,3)
                axs[binIndex].set_ylim(-2,3)
            
            plt.savefig('{directory}binnedBeamAlign_goTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! NOGO TONE !!
        if keyDataSheet['exp_seshPhase'] > 2:
            fig,axs=plt.subplots(nBins,1)
            axs[0].set_title('nogo tone - {}'.format(sites[s]))
            for binIndex in range(nBins):
                vline1=axs[binIndex].axvline(x=0, color=(1.0, 0.647, 0.0), linestyle='--', label='start cue')
                vline2=axs[binIndex].axvline(x=2.5, color='k', linestyle='--', label='tone')
                
                plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_nogoToneOnly'][s][binIndex], 'r', label='cue-align DA')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_DA_mean_nogoToneOnly'][s][binIndex] - alignmentOutputBinned['cueAlign_DA_std_nogoToneOnly'][s][binIndex], alignmentOutputBinned['cueAlign_DA_mean_nogoToneOnly'][s][binIndex] + alignmentOutputBinned['cueAlign_DA_std_nogoToneOnly'][s][binIndex], color='r', alpha=0.3)
                plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_nogoToneOnly'][s][binIndex], 'g', label='cue-align ACh')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['cueAlign_ACh_mean_nogoToneOnly'][s][binIndex] - alignmentOutputBinned['cueAlign_ACh_std_nogoToneOnly'][s][binIndex], alignmentOutputBinned['cueAlign_ACh_mean_nogoToneOnly'][s][binIndex] + alignmentOutputBinned['cueAlign_ACh_std_nogoToneOnly'][s][binIndex], color='g', alpha=0.3)
            
                axs[binIndex].set_xlabel('time from start (seconds)')
                axs[binIndex].set_ylabel('dF/F (%)')
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align DA'),
                         Line2D([binIndex], [binIndex], color='g', lw=2, label='cue-align ACh'),
                                Line2D([binIndex], [binIndex], color=(1.0, 0.647, 0.0), linestyle='--', lw=2, label='start cue'),
                                Line2D([binIndex], [binIndex], color='k', linestyle='--', lw=2, label='tone')]
                axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
                axs[binIndex].set_xlim(-2,6)
                axs[binIndex].set_ylim(-2,3)
        
            plt.savefig('{directory}binnedCueAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
            
            ### !! ALIGN TO BEAM 
            fig,axs=plt.subplots(nBins,1)
            axs[0].set_title('beam break after nogo Tone - {}'.format(sites[s]))
            for binIndex in range(nBins):
                vline1=axs[binIndex].axvline(x=0, color='b', linestyle='--', label='beam break')
            
                plot1=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean_nogoToneOnly'][s][binIndex], 'r', label='beam-align DA')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_DA_mean_nogoToneOnly'][s][binIndex] - alignmentOutputBinned['beamAlign_DA_std_nogoToneOnly'][s][binIndex], alignmentOutputBinned['beamAlign_DA_mean_nogoToneOnly'][s][binIndex] + alignmentOutputBinned['beamAlign_DA_std_nogoToneOnly'][s][binIndex], color='r', alpha=0.3)
                plot2=axs[binIndex].plot(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean_nogoToneOnly'][s][binIndex], 'g', label='beam-align ACh')
                axs[binIndex].fill_between(alignmentOutput['xAxisPoints'], alignmentOutputBinned['beamAlign_ACh_mean_nogoToneOnly'][s][binIndex] - alignmentOutputBinned['beamAlign_ACh_std_nogoToneOnly'][s][binIndex], alignmentOutputBinned['beamAlign_ACh_mean_nogoToneOnly'][s][binIndex] + alignmentOutputBinned['beamAlign_ACh_std_nogoToneOnly'][s][binIndex], color='g', alpha=0.3)
            
                axs[binIndex].set_xlabel('time from beam break (seconds)')
                axs[binIndex].set_ylabel('dF/F (%)')
                lines = [Line2D([binIndex], [binIndex], color='r', lw=2, label='beam-align DA'),
                         Line2D([binIndex], [binIndex], color='r', lw=2, label='cue-align ACh'),
                                Line2D([binIndex], [binIndex], color='b', linestyle='--', lw=2, label='beam break')]
                axs[binIndex].legend(handles=lines, loc='upper right', bbox_to_anchor=(0.95, 0.98))
                axs[binIndex].set_xlim(-2,3)
                axs[binIndex].set_ylim(-2,3)
        
            plt.savefig('{directory}binnedBeamAlign_nogoTone_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, hemi=sites[s], exp=keyDataSheet['experiment'], sesh=keyDataSheet['exp_seshNumber'], mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            #plt.close()
                    
    #%% Part 2.0 - return / save
    
    print(' Done.')
    # File path
    file_path = 'alignedDFOF_{}.pickle'.format(rxSite)
    # Writing the dictionary to a file using pickle
    with open(os.path.join(rxPath, file_path), 'wb') as file:
        pickle.dump({'alignmentOutput': alignmentOutput,'alignmentOutputBinned': alignmentOutputBinned,'alignedData': alignedData}, file)
    
    return alignmentOutput