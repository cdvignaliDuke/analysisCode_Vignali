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
import glob
import pickle
import tdt # tdt library
import sys
import glob
import pandas as pd
#

os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
import localFunctions as lf
import rxPreproc 
import rxDataOrg
import rxEpocAlign
#%% Part 0.1 - logical chain

userInputs = {'doCohort':[], 'doRxOnly':[], 'pullMPC':[], 'pullSynapse':[], 'errorPullMPC':[]}
userInputs['doCohort'] = bool(int(input('Do you want this script to iterate through every mice in this cohort?: ')))

if not userInputs['doCohort']:
# !! SET LOGIC CHAIN !!
    print('These will be logical decisions, so type either 1 or 0')
    userInputs['doRxOnly'] = True #bool(int(input('Do you want to analyze ONLY the recording data?: ')))
    userInputs['oneMouse'] = True #bool(int(input('Do you want to analyze ONE mouse [alt. is group analysis]?: '))) #    if not userInputs['doRxOnly']:
    ## !! THIS ASSUMES YOU HAVE ALREADY RAN BX DATA THROUGH bxParent !!
    userInputs['pullSynapse'] = True #bool(int(input('Do you have fiber photometry data associated with this MPC file?: ')))
    userInputs['pullMPC'] = False #bool(int(input('Do you need to preprocess the MPC file(s)?: '))) #    else: #        userInputs['pullMPC'] = True    
    
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
        
elif userInputs['doCohort']:
    userInputs['doRxOnly'] = True
    userInputs['oneMouse'] = False
    # !! SET PATH !!
    inputPath = ('S:\\Private\\Data\\Vignali cv105\\data\\MARS\\')
    outputPath = ('S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\')
    # !! SET MICE !!
    with open(inputPath + '\\' + 'activeCohort.txt','r') as dataLines:
        allMice=[]
        # !! ADD DATE !!
        for line in dataLines:
            allMice.append(line[0:3])
                
#%% Part 1.0 - for single mouse - format data folder and set root path

if userInputs['oneMouse']:
    # !! MOVE FILES !!
    lf.moveFiles(thisMouse,inputPath,suffix='-')
#    if date == []:
#       date = input('Date not found: input experiment date [YYMMDD]: ') 
    
    # !! SET DATE !!
    allFolders = os.listdir(os.path.join(inputPath,thisMouse)); currFolders = [fol for fol in allFolders if (fol.startswith('25') or fol.startswith('24') or fol.startswith('23'))] #REMOVED os.path.isdir(inputPath + fol + '\\') and
    print(currFolders)
    thisDay = input('Which day from currFolders: ') 
    
    if int(thisMouse) <= 18:
        rxSite = input('What hemisphere is recording data from [left / right / both]: ')
    else:
        rxSite = 'both'
        
    # !! SET PATH !!
    dataPath = inputPath + thisMouse + '\\' + thisDay + '\\'
    bxPath = outputPath + thisMouse + '\\bx\\'
    rxPath = outputPath + thisMouse + '\\' + thisDay + '\\fiber\\'
    os.makedirs(os.path.dirname(rxPath), exist_ok=True)

    toLoad = [file for file in os.listdir(dataPath)] #and os.path.isfile(os.path.join(inputPath, file))]
    for file in range(len(toLoad)):
        if os.path.isdir(os.path.join(dataPath,toLoad[file])):
            rxDataFolder = toLoad[file] + '\\'

    #%% Part 1.1 - import bx data
    
    if userInputs['pullMPC']:
        import bxParent
    elif not userInputs['pullMPC']:
        print('Importing bx data: ')
        import pickle
        bxDataSheet = pickle.load(open(os.path.join(bxPath, 'suppDataFile.pkl'), 'rb'))
    keyDataSheet = {}
    dayIndex = [i for i, item in enumerate(currFolders) if item == thisDay][0]-1; 
    if dayIndex == -1: dayIndex = 0
    # raw values (for session wide plot)
    keyDataSheet['trialStart'] = bxDataSheet['rawStart'][dayIndex] + 0.5 #500 msec delay between house light turning off and trial light turning on in MedPC design   
    keyDataSheet['trialEnd'] = bxDataSheet['rawEnd'][dayIndex]
    keyDataSheet['trialResponses'] = bxDataSheet['tBeamBreak'][dayIndex]
    keyDataSheet['trialReward'] = bxDataSheet['rawReward'][dayIndex]
    keyDataSheet['trialFirstResponses'] = bxDataSheet['rawBeamBreak_calc'][dayIndex]
    keyDataSheet['trailNotFirstResponses'] = bxDataSheet['rawBeamBreak_calcOthers'][dayIndex]
    
    # values aligned to trial start (for trial averaged plots)
    keyDataSheet['startAligned_firstResponse'] = bxDataSheet['tBeamBreak_calc'][dayIndex]
    keyDataSheet['startAligned_reward'] = bxDataSheet['tReward'][dayIndex]
    keyDataSheet['startAligned_notFirstresponses'] = bxDataSheet['tBeamBreak_calcOthers'][dayIndex]
    
    # trial indices
    keyDataSheet['trialType'] = bxDataSheet['tType'][dayIndex]
    keyDataSheet['accOfResponse'] = bxDataSheet['responseMatrix'][dayIndex]
    keyDataSheet['trialType_calc'] = bxDataSheet['tType_calc'][dayIndex]
    
    #ensure calc dict is full length
    if len(keyDataSheet['trialType_calc']) < len(keyDataSheet['trialType']):
        keyDataSheet['trialType_calc'] = np.append(keyDataSheet['trialType_calc'], np.zeros(len(keyDataSheet['trialType']) - len(keyDataSheet['trialType_calc'])))
    
    # save info
    keyDataSheet['experiment'] = bxDataSheet['allDataSheet'].at[dayIndex,'Experiment'] #fiber during what experiment?
    keyDataSheet['sessionNumber'] = bxDataSheet['allDataSheet'].loc[bxDataSheet['allDataSheet']['Date'] == thisDay].index[0]
    keyDataSheet['exp_seshNumber'] = bxDataSheet['phaseIndex'][:keyDataSheet['sessionNumber'] + 1].count(bxDataSheet['phaseIndex'][keyDataSheet['sessionNumber']])
    keyDataSheet['exp_seshPhase'] = bxDataSheet['phaseIndex'][keyDataSheet['sessionNumber']]
    keyDataSheet['thisMouse'] = thisMouse #in what mouse
    keyDataSheet['thisDay'] = thisDay #on which day
    keyDataSheet['rxSite'] = rxSite #from what region (as of 05/13/2024: left = DLS, right = DMS)
                
    #%% Part 1.2 - fiber data organization and preprocessing
    
    fibDataSaved, epocSaved = rxDataOrg.loadRawData(userInputs, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse, thisDay, currFolders, rxSite)
    
    for key, value in epocSaved.items():
        globals()[key] = value
    for key, value in fibDataSaved.items():
            globals()[key] = value
            
    #%% Part 1.3 - preprocessing section [adapted from 'photometry data preprocessing.ipynb' in sourceCode folder]
    pickles = [file for file in os.listdir(rxPath) if file.endswith('.pickle')]
    sites = ['left','right']
    if pickles == []:
        hemisPreproc = {} # dictionary for consolidation of data
        if rxSite == 'left':
            preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
            hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
        elif rxSite == 'right':
            preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
            hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
        elif rxSite == 'both':
            preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'left', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
            hemisPreproc = rxPreproc.makeHemisDir('left', hemisPreproc, preprocOutput, rxPath)
            preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'right', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
            hemisPreproc = rxPreproc.makeHemisDir('right', hemisPreproc, preprocOutput, rxPath)
    else:    
        if 'hemisPreproc_{}.pickle'.format(rxSite) in pickles:
            if rxSite == 'both':
                hemisPreproc = {} # dictionary for consolidation of data
                for i in range(len(sites)):
                    file_path = 'hemisPreproc_{}.pickle'.format(sites[i])
                    if file_path in pickles:
                        if os.path.getsize(os.path.join(rxPath, file_path))>0:
                            print('   Saved preprocessing files, pulling now: ')
                            with open(os.path.join(rxPath, file_path), 'rb') as file:
                                hemisPreproc = pickle.load(file)
                            #hemisPreproc = rxPreproc.makeHemisDir(sites[i], hemisPreproc, preprocOutput, rxPath)
            else:
                hemisPreproc = {} # dictionary for consolidation of data
                file_path = 'hemisPreproc_{}.pickle'.format(rxSite)
                if file_path in pickles:
                    if os.path.getsize(os.path.join(rxPath, file_path))>0:
                        print('   Saved preprocessing files, pulling now: ')
                        with open(os.path.join(rxPath, file_path), 'rb') as file:
                            hemisPreproc = pickle.load(file)
                        #hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
        else:
            hemisPreproc = {} # dictionary for consolidation of data
            if rxSite == 'left':
                preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
            elif rxSite == 'right':
                preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
            elif rxSite == 'both':
                preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'left', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir('left', hemisPreproc, preprocOutput, rxPath)
                preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'right', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir('right', hemisPreproc, preprocOutput, rxPath)
    #%% Part 1.4 - cue and beam alignment
    if rxSite == 'left':
        alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
    elif rxSite == 'right':
        alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
    elif rxSite == 'both':
        alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
    
    #%% Part 1.4.1 - pull raw fluorescence data
    folders = os.listdir(outputPath + thisMouse) #from path, list folders
    currFolders = folders[:-1]
    print(currFolders)
    
    laserDataSheet = pd.DataFrame()
    laserDataSheet = pd.DataFrame(columns = ['left_isosRaw_first500', 'left_gAChRaw_first500', 'left_rdLightRaw_first500', 'left_isosRaw_last500', 'left_gAChRaw_last500', 'left_rdLightRaw_last500', 'right_isosRaw_first500', 'right_gAChRaw_first500', 'right_rdLightRaw_first500', 'right_isosRaw_last500', 'right_gAChRaw_last500', 'right_rdLightRaw_last500', 'left_isosRaw_first500_mean', 'left_gAChRaw_first500_mean', 'left_rdLightRaw_first500_mean', 'left_isosRaw_last500_mean', 'left_gAChRaw_last500_mean', 'left_rdLightRaw_last500_mean', 'right_isosRaw_first500_mean', 'right_gAChRaw_first500_mean', 'right_rdLightRaw_first500_mean', 'right_isosRaw_last500_mean', 'right_gAChRaw_last500_mean', 'right_rdLightRaw_last500_mean', 'left_isosRaw_first500_std', 'left_gAChRaw_first500_std', 'left_rdLightRaw_first500_std', 'left_isosRaw_last500_std', 'left_gAChRaw_last500_std', 'left_rdLightRaw_last500_std', 'right_isosRaw_first500_std', 'right_gAChRaw_first500_std', 'right_rdLightRaw_first500_std', 'right_isosRaw_last500_std', 'right_gAChRaw_last500_std', 'right_rdLightRaw_last500_std'])
    
    for fol in currFolders:
        # Reading the dictionary from a file using pickle
        musHemisPreproc = {}; 
        for s in range(len(glob.glob(os.path.join(outputPath + thisMouse + '\\' + fol + '\\fiber\\', 'hemisPreproc_*.pickle')))):
            loaded = pickle.load(open(glob.glob(os.path.join(outputPath + thisMouse + '\\' + fol + '\\fiber\\', 'hemisPreproc_*.pickle'))[s], 'rb'))
            if int(thisMouse) <= 18:
                musHemisPreproc = loaded
            else:
                musHemisPreproc[sites[s]] = loaded[sites[s]]

        if 'left' in musHemisPreproc:
            laserDataSheet.at[fol,'left_isosRaw_first500'] = (musHemisPreproc['left']['isos_raw'][200:700])
            laserDataSheet.at[fol,'left_gAChRaw_first500'] = (musHemisPreproc['left']['gACh_raw'][200:700])
            laserDataSheet.at[fol,'left_rdLightRaw_first500'] = (musHemisPreproc['left']['rdLight_raw'][200:700])
            laserDataSheet.at[fol,'left_isosRaw_last500'] = (musHemisPreproc['left']['isos_raw'][-500:])
            laserDataSheet.at[fol,'left_gAChRaw_last500'] = (musHemisPreproc['left']['gACh_raw'][-500:])
            laserDataSheet.at[fol,'left_rdLightRaw_last500'] = (musHemisPreproc['left']['rdLight_raw'][-500:])
            laserDataSheet.at[fol,'left_isosRaw_first500_mean'] = np.mean(musHemisPreproc['left']['isos_raw'][200:700]); laserDataSheet.at[fol,'left_isosRaw_first500_std'] = np.std(musHemisPreproc['left']['isos_raw'][200:700])
            laserDataSheet.at[fol,'left_gAChRaw_first500_mean'] = np.mean(musHemisPreproc['left']['gACh_raw'][200:700]); laserDataSheet.at[fol,'left_gAChRaw_first500_std'] = np.std(musHemisPreproc['left']['gACh_raw'][200:700])
            laserDataSheet.at[fol,'left_rdLightRaw_first500_mean'] = np.mean(musHemisPreproc['left']['rdLight_raw'][200:700]); laserDataSheet.at[fol,'left_rdLightRaw_first500_std'] = np.std(musHemisPreproc['left']['rdLight_raw'][200:700])
            laserDataSheet.at[fol,'left_isosRaw_last500_mean'] = np.mean(musHemisPreproc['left']['isos_raw'][-500:]); laserDataSheet.at[fol,'left_isosRaw_last500_std'] = np.std(musHemisPreproc['left']['isos_raw'][-500:])
            laserDataSheet.at[fol,'left_gAChRaw_last500_mean'] = np.mean(musHemisPreproc['left']['gACh_raw'][-500:]); laserDataSheet.at[fol,'left_gAChRaw_last500_std'] = np.std(musHemisPreproc['left']['gACh_raw'][-500:])
            laserDataSheet.at[fol,'left_rdLightRaw_last500_mean'] = np.mean(musHemisPreproc['left']['rdLight_raw'][-500:]); laserDataSheet.at[fol,'left_rdLightRaw_last500_std'] = np.std(musHemisPreproc['left']['rdLight_raw'][-500:])
        if 'right' in musHemisPreproc:
            laserDataSheet.at[fol,'right_isosRaw_first500'] = (musHemisPreproc['right']['isos_raw'][200:700])
            laserDataSheet.at[fol,'right_gAChRaw_first500'] = (musHemisPreproc['right']['gACh_raw'][200:700])
            laserDataSheet.at[fol,'right_rdLightRaw_first500'] = (musHemisPreproc['right']['rdLight_raw'][200:700])
            laserDataSheet.at[fol,'right_isosRaw_last500'] = (musHemisPreproc['right']['isos_raw'][-500:])
            laserDataSheet.at[fol,'right_gAChRaw_last500'] = (musHemisPreproc['right']['gACh_raw'][-500:])
            laserDataSheet.at[fol,'right_rdLightRaw_last500'] = (musHemisPreproc['right']['rdLight_raw'][-500:])
            laserDataSheet.at[fol,'right_isosRaw_first500_mean'] = np.mean(musHemisPreproc['right']['isos_raw'][200:700]); laserDataSheet.at[fol,'right_isosRaw_first500_std'] = np.std(musHemisPreproc['right']['isos_raw'][200:700])
            laserDataSheet.at[fol,'right_gAChRaw_first500_mean'] = np.mean(musHemisPreproc['right']['gACh_raw'][200:700]); laserDataSheet.at[fol,'right_gAChRaw_first500_std'] = np.std(musHemisPreproc['right']['gACh_raw'][200:700])
            laserDataSheet.at[fol,'right_rdLightRaw_first500_mean'] = np.mean(musHemisPreproc['right']['rdLight_raw'][200:700]); laserDataSheet.at[fol,'right_rdLightRaw_first500_std'] = np.std(musHemisPreproc['right']['rdLight_raw'][200:700])
            laserDataSheet.at[fol,'right_isosRaw_last500_mean'] = np.mean(musHemisPreproc['right']['isos_raw'][-500:]); laserDataSheet.at[fol,'right_isosRaw_last500_std'] = np.std(musHemisPreproc['right']['isos_raw'][-500:])
            laserDataSheet.at[fol,'right_gAChRaw_last500_mean'] = np.mean(musHemisPreproc['right']['gACh_raw'][-500:]); laserDataSheet.at[fol,'right_gAChRaw_last500_std'] = np.std(musHemisPreproc['right']['gACh_raw'][-500:])
            laserDataSheet.at[fol,'right_rdLightRaw_last500_mean'] = np.mean(musHemisPreproc['right']['rdLight_raw'][-500:]); laserDataSheet.at[fol,'right_rdLightRaw_last500_std'] = np.std(musHemisPreproc['right']['rdLight_raw'][-500:])
        
#      
    plt.figure(figsize=(20, 16))
    plt.title('Raw fluorescence at start/end of session - DLS',fontsize=36); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Power (V)',fontsize=36)
    plt.tick_params(axis='both', labelsize=48); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist(), color='blue', linewidth=10, markersize=15, marker='o', label='Isos_start')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['left_isosRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='blue')
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_gAChRaw_first500_mean'].dropna().tolist(), color='green', linewidth=10, markersize=15, marker='o', label='gACh_start')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_gAChRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['left_gAChRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='green')
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_rdLightRaw_first500_mean'].dropna().tolist(), color='red', linewidth=10, markersize=15, marker='o', label='rdLi_start')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['left_rdLightRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['left_rdLightRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='red')
#
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist(), color='blue', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='Isos_end')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['left_isosRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='blue')
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_gAChRaw_last500_mean'].dropna().tolist(), color='green', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='GACh_end')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_gAChRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['left_gAChRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='green')
    plt.plot((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_rdLightRaw_last500_mean'].dropna().tolist(), color='red', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='RdLi_end')
    plt.errorbar((range(1, len(laserDataSheet['left_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['left_rdLightRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['left_rdLightRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='red')

    plt.legend(fontsize=18, title='session', loc='lower right', bbox_to_anchor=(1,.05));
    plt.savefig('{directory}rawFluor_start-end_leftHemi.pdf'.format(directory=rxPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); #plt.close()

#~~
    plt.figure(figsize=(20, 16))
    plt.title('Raw fluorescence at start/end of session - DMS',fontsize=36); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Power (V)',fontsize=36)
    plt.tick_params(axis='both', labelsize=48); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist(), color='blue', linewidth=10, markersize=15, marker='o', label='Isos_start')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['right_isosRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='blue')
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_gAChRaw_first500_mean'].dropna().tolist(), color='green', linewidth=10, markersize=15, marker='o', label='gACh_start')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_gAChRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['right_gAChRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='green')
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_rdLightRaw_first500_mean'].dropna().tolist(), color='red', linewidth=10, markersize=15, marker='o', label='rDLi_start')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_first500_mean'].dropna().tolist())+1)), laserDataSheet['right_rdLightRaw_first500_mean'].dropna().tolist(), yerr=laserDataSheet['right_rdLightRaw_first500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='red')
#
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist(), color='blue', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='Isos_end')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['right_isosRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='blue')
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_gAChRaw_last500_mean'].dropna().tolist(), color='green', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='gACh_end')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_gAChRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['right_gAChRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='green')
    plt.plot((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_rdLightRaw_last500_mean'].dropna().tolist(), color='red', linewidth=10, markersize=15, marker='o', linestyle='dashed', label='rDLi_end')
    plt.errorbar((range(1, len(laserDataSheet['right_isosRaw_last500_mean'].dropna().tolist())+1)), laserDataSheet['right_rdLightRaw_last500_mean'].dropna().tolist(), yerr=laserDataSheet['right_rdLightRaw_last500_std'].dropna().tolist(), fmt='none', capsize=7.5, elinewidth=3, color='red')
    
    plt.legend(fontsize=18, title='session', loc='lower right', bbox_to_anchor=(1,.75));
    plt.savefig('{directory}rawFluor_start-end_rightHemi.pdf'.format(directory=rxPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); #plt.close()
          
                  
#%% Part 2.0 - for all mice - format data folder and set root path

elif not userInputs['oneMouse']:
    
    thisDay = (input('Which day from currFolders: '))
    
    rxSite = input('What hemisphere is recording data from [left / right / both]: ')
    
    for mus in range(len(allMice)):
        thisMouse = allMice[mus]
        plt.close('all'); print('Starting mouse {}:'.format(thisMouse))


        allFolders = os.listdir(os.path.join(inputPath,thisMouse)); currFolders = [fol for fol in allFolders if (fol.startswith('25') or fol.startswith('24') or fol.startswith('23'))] #REMOVED os.path.isdir(inputPath + fol + '\\') and
        # !! MOVE FILES !!
        lf.moveFiles(thisMouse,inputPath,suffix='-')
    
        # !! SET PATH !!
        dataPath = inputPath + thisMouse + '\\' + thisDay + '\\'
        bxPath = outputPath + thisMouse + '\\bx\\'
        rxPath = outputPath + thisMouse + '\\' + thisDay + '\\fiber\\'
        os.makedirs(os.path.dirname(rxPath), exist_ok=True)
    
        toLoad = [file for file in os.listdir(dataPath)] #and os.path.isfile(os.path.join(inputPath, file))]
        for file in range(len(toLoad)):
            if os.path.isdir(os.path.join(dataPath,toLoad[file])):
                rxDataFolder = toLoad[file]
            
        #%% Part 2.1 - import bx data
        
        if userInputs['pullMPC']:
            import bxParent
        elif not userInputs['pullMPC']:
            print('Importing bx data: ')
            import pickle
            bxDataSheet = pickle.load(open(os.path.join(bxPath, 'suppDataFile.pkl'), 'rb'))
        keyDataSheet = {}
        dayIndex = [i for i, item in enumerate(currFolders) if item == thisDay][0]-1
        # raw values (for session wide plot)
        keyDataSheet['trialStart'] = bxDataSheet['rawStart'][dayIndex] + 0.5 #500 msec delay between house light turning off and trial light turning on in MedPC design   
        keyDataSheet['trialEnd'] = bxDataSheet['rawEnd'][dayIndex]
        keyDataSheet['trialResponses'] = bxDataSheet['tBeamBreak'][dayIndex]
        keyDataSheet['trialReward'] = bxDataSheet['rawReward'][dayIndex]
        keyDataSheet['trialFirstResponses'] = bxDataSheet['rawBeamBreak_calc'][dayIndex]
        keyDataSheet['trailNotFirstResponses'] = bxDataSheet['rawBeamBreak_calcOthers'][dayIndex]
        
        # values aligned to trial start (for trial averaged plots)
        keyDataSheet['startAligned_firstResponse'] = bxDataSheet['tBeamBreak_calc'][dayIndex]
        keyDataSheet['startAligned_reward'] = bxDataSheet['tReward'][dayIndex]
        keyDataSheet['startAligned_notFirstresponses'] = bxDataSheet['tBeamBreak_calcOthers'][dayIndex]
        
        # trial indices
        keyDataSheet['trialType'] = bxDataSheet['tType'][dayIndex]
        keyDataSheet['accOfResponse'] = bxDataSheet['responseMatrix'][dayIndex]
        keyDataSheet['trialType_calc'] = bxDataSheet['tType_calc'][dayIndex]
        
        #ensure calc dict is full length
        if len(keyDataSheet['trialType_calc']) < len(keyDataSheet['trialType']):
            keyDataSheet['trialType_calc'] = np.append(keyDataSheet['trialType_calc'], np.zeros(len(keyDataSheet['trialType']) - len(keyDataSheet['trialType_calc'])))
        
        # save info
        keyDataSheet['experiment'] = bxDataSheet['allDataSheet'].at[dayIndex,'Experiment'] #fiber during what experiment?
        keyDataSheet['sessionNumber'] = bxDataSheet['allDataSheet'].loc[bxDataSheet['allDataSheet']['Date'] == thisDay].index[0]
        keyDataSheet['exp_seshNumber'] = bxDataSheet['phaseIndex'][:keyDataSheet['sessionNumber'] + 1].count(bxDataSheet['phaseIndex'][keyDataSheet['sessionNumber']])
        keyDataSheet['exp_seshPhase'] = bxDataSheet['phaseIndex'][keyDataSheet['sessionNumber']]
        keyDataSheet['thisMouse'] = thisMouse #in what mouse
        keyDataSheet['thisDay'] = thisDay #on which day
        keyDataSheet['rxSite'] = rxSite #from what region (as of 05/13/2024: left = DLS, right = DMS)
                    
        #%% Part 2.2 - fiber data organization and preprocessing
        
        fibDataSaved, epocSaved = rxDataOrg.loadRawData(userInputs, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse, thisDay, currFolders, rxSite)
        
        for key, value in epocSaved.items():
            globals()[key] = value
        for key, value in fibDataSaved.items():
                globals()[key] = value
                
        #%% Part 2.3 - preprocessing section [adapted from 'photometry data preprocessing.ipynb' in sourceCode folder]
        pickles = [file for file in os.listdir(rxPath) if file.endswith('.pickle')]
        sites = ['left','right']
        if pickles == []:
            hemisPreproc = {} # dictionary for consolidation of data
            if rxSite == 'left':
                preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
            elif rxSite == 'right':
                preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
            elif rxSite == 'both':
                preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'left', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir('left', hemisPreproc, preprocOutput, rxPath)
                preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'right', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                hemisPreproc = rxPreproc.makeHemisDir('right', hemisPreproc, preprocOutput, rxPath)
        else:    
            if 'hemisPreproc_{}.pickle'.format(rxSite) in pickles:
                if rxSite == 'both':
                    hemisPreproc = {} # dictionary for consolidation of data
                    for i in range(len(sites)):
                        file_path = 'hemisPreproc_{}.pickle'.format(sites[i])
                        if file_path in pickles:
                            if os.path.getsize(os.path.join(rxPath, file_path))>0:
                                print('   Saved preprocessing files, pulling now: ')
                                with open(os.path.join(rxPath, file_path), 'rb') as file:
                                    hemisPreproc = pickle.load(file)
                                #hemisPreproc = rxPreproc.makeHemisDir(sites[i], hemisPreproc, preprocOutput, rxPath)
                else:
                    hemisPreproc = {} # dictionary for consolidation of data
                    file_path = 'hemisPreproc_{}.pickle'.format(rxSite)
                    if file_path in pickles:
                        if os.path.getsize(os.path.join(rxPath, file_path))>0:
                            print('   Saved preprocessing files, pulling now: ')
                            with open(os.path.join(rxPath, file_path), 'rb') as file:
                                hemisPreproc = pickle.load(file)
                            #hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
            else:
                hemisPreproc = {} # dictionary for consolidation of data
                if rxSite == 'left':
                    preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                    hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
                elif rxSite == 'right':
                    preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, rxSite, keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                    hemisPreproc = rxPreproc.makeHemisDir(rxSite, hemisPreproc, preprocOutput, rxPath)
                elif rxSite == 'both':
                    preprocOutput = rxPreproc.preprocessing(data405_LH, data465_LH, data560_LH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'left', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                    hemisPreproc = rxPreproc.makeHemisDir('left', hemisPreproc, preprocOutput, rxPath)
                    preprocOutput = rxPreproc.preprocessing(data405_RH, data465_RH, data560_RH, fs, epoc_ticker, epoc_beamBreak, rxPath, 'right', keyDataSheet['experiment'], keyDataSheet['exp_seshNumber'], thisMouse)
                    hemisPreproc = rxPreproc.makeHemisDir('right', hemisPreproc, preprocOutput, rxPath)
                        
        #%% Part 2.4 - cue and beam alignment
        if rxSite == 'left':
            alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
        elif rxSite == 'right':
            alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
        elif rxSite == 'both':
            alignmentOutput = rxEpocAlign.alignment(hemisPreproc, epoc_onsetITI, epoc_onsetCues, keyDataSheet, epoc_beamBreak, rxPath, rxSite, thisMouse, bxDataSheet, dayIndex)
              
     
    
#%%
sys.exit()

#%%



    #Create a mean signal, standard error of signal, and DC offset
if rxSite == 'left':
    mean560_LH = np.mean(hemisPreproc['left']['rdLight_dF_F'], axis=0)
    std560_LH = np.std(hemisPreproc['left']['rdLight_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['left']['rdLight_dF_F']))
    dc560_LH = np.mean(mean560_LH)
    mean465_LH = np.mean(hemisPreproc['left']['gACh_dF_F'], axis=0)
    std465_LH = np.std(hemisPreproc['left']['gACh_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['left']['gACh_dF_F']))
    dc465_LH = np.mean(mean465_LH)
if rxSite == 'right':
    mean560_RH = np.mean(hemisPreproc['right']['rdLight_dF_F'], axis=0)
    std560_RH = np.std(hemisPreproc['right']['rdLight_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['right']['rdLight_dF_F']))
    dc560_RH = np.mean(mean560_RH)
    mean465_RH = np.mean(hemisPreproc['right']['gACh_dF_F'], axis=0)
    std465_RH = np.std(hemisPreproc['right']['gACh_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['right']['gACh_dF_F']))
    dc465_RH = np.mean(mean465_RH)
if rxSite == 'both':
    mean560_LH = np.mean(hemisPreproc['left']['rdLight_dF_F'], axis=0)
    std560_LH = np.std(hemisPreproc['left']['rdLight_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['left']['rdLight_dF_F']))
    dc560_LH = np.mean(mean560_LH)
    mean465_LH = np.mean(hemisPreproc['left']['gACh_dF_F'], axis=0)
    std465_LH = np.std(hemisPreproc['left']['gACh_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['left']['gACh_dF_F']))
    dc465_LH = np.mean(mean465_LH)
    mean560_RH = np.mean(hemisPreproc['right']['rdLight_dF_F'], axis=0)
    std560_RH = np.std(hemisPreproc['right']['rdLight_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['right']['rdLight_dF_F']))
    dc560_RH = np.mean(mean560_RH)
    mean465_RH = np.mean(hemisPreproc['right']['gACh_dF_F'], axis=0)
    std465_RH = np.std(hemisPreproc['right']['gACh_dF_F'], axis=0)/np.sqrt(len(hemisPreproc['right']['gACh_dF_F']))
    dc465_RH = np.mean(mean465_RH)


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
    
else:
    sys.exit()