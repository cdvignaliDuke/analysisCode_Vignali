# -*- coding: utf-8 -*-
"""
Created on Fri Dec 08 14:35:25 2023

@author: Carlo Vignali

@def: This code will contain 4 primary parts:
    1) read in MPC file and separate variables from datasheet into arrays
    2) when prompted, import data from DeepLabCut with timestamps for specific events
    3) import preprocessed fiber photometry data (or initialize Synapse to preprocess data using MPC timestamps)
    4) run primary analyses:
        this will have a number of parts that will follow when I actually have data to manipulate and consider
"""
#%%
# import local functions and open variables
import os
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
import localFunctions as lf

#%% Part 0 - determine subsections and early decision forks

if not 'userInputs' in locals():
    print('These will be logical decisions, so type either 1 or 0')
    userInputs = {'doBxOnly':[], 'pullMPC':[], 'pullSynapse':[], 'errorPullMPC':[]}
    userInputs['doBxOnly'] = bool(int(input('Do you want to analyze ONLY the behavioral data?: ')))
    userInputs['oneMouse'] = bool(int(input('Do you want to analyze one mouse or do group analysis?: ')))
    if not userInputs['doBxOnly']:
        userInputs['pullSynapse'] = bool(int(input('Do you have fiber photometry data associated with this MPC file?: ')))
        userInputs['pullMPC'] = bool(int(input('Do you need to preprocess the MPC file(s)?: ')))
    else:
        userInputs['pullMPC'] = True    

    inputPath = ('S:\\Private\\Data\\Vignali cv105\\data\\MARS\\')
    outputPath = ('S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\')
    
    if userInputs['oneMouse']:
        # !! SET MOUSE !!
        #determine path to datafile - experiment file + animal number + date
        #from path, list folders, then user provides input to choose the animal of interest
        currFolders = os.listdir(inputPath)[1:]
        print(currFolders)
        thisMouse = currFolders[int(input('Which mouse from currFolders: '))-1] # THIS LINE SUBTRACTS 1 FROM THE RESPONSE TO CORRECT INPUT, ASSUMING THE USER IS INPUTTING THE MOUSE ID OF INTEREST AND NOT THE PYTHON INDEX FOR THAT MOUSE
        
#%% Part 1 - read MPC
"""
Syntax:
    A- trial start
    B- trial end
    C- trial type (go = 1, no-go = 2, tooEarly (before tone initalized) = 0)
    D- trial ITI duration
    E- response window close
    F- beam break
    G- reward time
    H- 
    I- response counter (hit, correct rejection, miss, false alarm, too early, beam break)
    J- time played (should sync with C)
    K- trial result (trial-by-trial basis; should match up with other calculations based on beam break, reward and response window)
    L- 
    M- 
    N- max trials
    O- 
    P- 
    Q- 
    R- 
    S- 
    T- (0) session clock (end)
    U- 
    V- 
    W- trial (end)
    X- 
    Y- ITI range
    Z- 
"""
if userInputs['oneMouse']:
        
    # !! MOVE FILES !!
    lf.moveFiles(thisMouse,inputPath,fileType='.txt')
    
    #from path, list folders, then user provides input to choose the day of interest
    inputPath = inputPath + thisMouse + '\\'
    outputPath = outputPath + thisMouse + '\\bx\\'
    #ensure outputPath exists in the directory
    os.makedirs(os.path.dirname(outputPath), exist_ok=True)
    #
    if not userInputs['doBxOnly']:
        allFolders = os.listdir(inputPath); currFolders = [fol for fol in allFolders if os.path.isdir(os.path.join(inputPath, fol)) and (fol.startswith('24') or fol.startswith('23'))]
        print(currFolders)
        thisDay = currFolders[int(input('Which day from currFolders: '))]
    #final path for specific session analysis
        inputPath = inputPath + thisDay + '\\'
    
    # !! read MPC files into dataSheet !!
    # if all behavioral data is already preprocessed (MPC to excel sheets)
    allDataSheet = pd.DataFrame() #columns = ['Subject', 'Experiment', 'Date', 'startTime', 'endTime', 'totalTime', 'ITIRange', 'maxTrials', 'endTrial', 'trialStartTime', 'trialEndTime', 'trialType', 'trialITI', 'responseTime', 'rewardTime', 'responseCounter', 'toneTime', 'terminalBeamBreakTime', 'beamBreakTime', 'trialOutcome'])
    if not userInputs['pullMPC']:
        # Are we only going to look at behavioral data
        if userInputs['doBxOnly']:    
            allFolders = os.listdir(inputPath); currFolders = [fol for fol in allFolders if os.path.isdir(os.path.join(inputPath, fol)) and (fol.startswith('24') or fol.startswith('23'))]
            print(currFolders)
            dayDataArray=[]
            # if so, look though each folder for that mouse
            for day in currFolders:
                # if the folder is empty, tag it and move on
                if not os.listdir(inputPath + currFolders[day] +'\\'):
                    dayDataArray.append(False)
                else:
                    dayDataArray.append(True)
                    excelExtensions = ['.xlsx','.xls','.xlsm']
                    # look to see if there is an excel sheet for that day
                    excelFiles = [file for file in os.listdir(inputPath + currFolders[day]) if file.endswith(tuple(excelExtensions))]
                    if not excelFiles == []:
                        dataSheet = pd.read_excel(inputPath + currFolders[day])
                    # if there is not an excel sheet, return error and prompt user to either pause to convert MPC>excel or close script
                    else:
                        userInputs['errorPullMPC']=bool(input('The MPC data for this day is unprocessed, would you like to pause and export the MPC to an excel?: '))        
                        if userInputs['errorPullMPC']:
                            txtFiles = [file for file in os.listdir(inputPath + currFolders[day]) if file.endswith(tuple('.txt'))]
                            if txtFiles:
                                dataSheet = lf.exportMPC(inputPath + currFolders[day])
                            else:
                                print('Session - {} - will be skipped. There is no txt file in this folder.'.format[day])    
                        else:
                            print('Closing script. Export MPC file before restarting this script.')
                            sys.exit()
                allDataSheet = pd.concat([allDataSheet, dataSheet], ignore_index=True)            
        else: 
            dataSheet = pd.read_excel(inputPath)
            allDataSheet = pd.concat([allDataSheet, dataSheet], ignore_index=True)
    else:
        if userInputs['doBxOnly']:    
            allFolders = os.listdir(inputPath); currFolders = [fol for fol in allFolders if os.path.isdir(os.path.join(inputPath, fol)) and (fol.startswith('24') or fol.startswith('23'))]
            print(currFolders); 
            dayDataArray=[]
            # if so, look though each folder for that mouse
            for day in range(len(currFolders)):
                print(day); dataPath = inputPath + currFolders[day]
                # if the folder is empty, tag it and move on
                if not os.listdir(dataPath + '\\'):
                    dayDataArray.append(False)
                else: 
                    dayDataArray.append(True)
                    dataSheet = lf.exportMPC(dataPath)
                    allDataSheet = pd.concat([allDataSheet, dataSheet], ignore_index=True)
    
    
#%% Part 2 - analyze data for each day    
    
    rateAccuracy = [[] for _ in range(len(allDataSheet))]; rateAccuracyGo = [[] for _ in range(len(allDataSheet))]; rateAccuracyNo = [[] for _ in range(len(allDataSheet))]; rateTooEarly = [[] for _ in range(len(allDataSheet))];
    nTrials = [[] for _ in range(len(allDataSheet))]; maxTrials = [[] for _ in range(len(allDataSheet))]; rawStart = [[] for _ in range(len(allDataSheet))]; rawEnd = [[] for _ in range(len(allDataSheet))]; rawResponse = [[] for _ in range(len(allDataSheet))]; rawTone = [[] for _ in range(len(allDataSheet))]; rawReward = [[] for _ in range(len(allDataSheet))];
    frameRate = [[] for _ in range(len(allDataSheet))]; tStart = [[] for _ in range(len(allDataSheet))]; tEnd = [[] for _ in range(len(allDataSheet))]; tResponse = [[] for _ in range(len(allDataSheet))]; tTone = [[] for _ in range(len(allDataSheet))]; tReward = [[] for _ in range(len(allDataSheet))]; tTerminalBeamBreak = [[] for _ in range(len(allDataSheet))]; tBeamBreak = [[] for _ in range(len(allDataSheet))];
    responseCounter = [[] for _ in range(len(allDataSheet))]; tType = [[] for _ in range(len(allDataSheet))]; tType_calc = [[] for _ in range(len(allDataSheet))];
    responseMatrix = [[] for _ in range(len(allDataSheet))]; hits = [[] for _ in range(len(allDataSheet))]; nHit = [[] for _ in range(len(allDataSheet))]; misses = [[] for _ in range(len(allDataSheet))]; nMiss = [[] for _ in range(len(allDataSheet))]; correctRejects = [[] for _ in range(len(allDataSheet))]; nCorrectReject = [[] for _ in range(len(allDataSheet))]; falseAlarms = [[] for _ in range(len(allDataSheet))]; nFalseAlarm = [[] for _ in range(len(allDataSheet))]; tooEarlies = [[] for _ in range(len(allDataSheet))]; nTooEarly = [[] for _ in range(len(allDataSheet))]; missTrials = [[] for _ in range(len(allDataSheet))]; nMissTrial = [[] for _ in range(len(allDataSheet))]; phaseIndex = [[] for _ in range(len(allDataSheet))];
        
    for day in range(len(allDataSheet)):
        nTrials[day]=(int(allDataSheet.loc[day,'endTrial']))
        maxTrials[day]=(int(allDataSheet.loc[day,'maxTrials']))
        
        rawStart[day]=((allDataSheet.loc[day,'trialStartTime']))
        rawEnd[day]=((allDataSheet.loc[day,'trialEndTime']))
        # quick modification that adds the end session time to the last trial of rawEnd if the session reaches 3600 sec before the last trial ends
        if 0 in rawEnd[day]:
            #lastStart=min(min(np.where(rawStart[day] == 0)))-1
            #if rawEnd[day][lastStart] == 0:
            #    rawEnd[day][lastStart] = 3600.00
            rawEnd[day][np.where(rawEnd[day]==0)]=3600.00            
            
        rawResponse[day]=((allDataSheet.loc[day,'responseTime']))
        rawTone[day]=((allDataSheet.loc[day,'toneTime']))
        rawReward[day]=((allDataSheet.loc[day,'rewardTime']))
        
        frameRate[day]=(30)
        tStart[day]=(np.zeros(nTrials[day]))
        for trial in range(len(tStart[day])):
            maths = rawEnd[day][trial] - rawStart[day][trial]  
            tEnd[day] = np.append(tEnd[day], maths)
        for trial in range(len(tStart[day])):
            maths = rawResponse[day][trial] - rawStart[day][trial]  
            tResponse[day] = np.append(tResponse[day], maths)
        for trial in range(len(tStart[day])):
            maths = rawReward[day][trial] - rawStart[day][trial]  
            tReward[day] = np.append(tReward[day], maths)
        if np.count_nonzero(rawTone[day]) > 0:    
            for trial in range(len(tStart[day])):
                maths = rawTone[day][trial] - rawStart[day][trial]  
                tTone[day] = np.append(tTone[day], maths)
        else:
            tTone[day] = rawTone[day]
        
        responseCounter[day] = allDataSheet.loc[day,'responseCounter']   
        responseMatrix[day] = allDataSheet.loc[day,'trialOutcome']
        if isinstance(responseMatrix[day], np.ndarray):
            hits[day] = [1 if x == 1 else 0 for x in responseMatrix[day]]
            nHit[day] = sum(hits[day])
            misses[day] = [1 if x == 3 else 0 for x in responseMatrix[day]]
            nMiss[day] = sum(misses[day])
            correctRejects[day] = [1 if x == 2 else 0 for x in responseMatrix[day]]
            nCorrectReject[day] = sum(correctRejects[day])
            falseAlarms[day] = [1 if x == 4 else 0 for x in responseMatrix[day]]
            nFalseAlarm[day] = sum(falseAlarms[day])
            tooEarlies[day] = [1 if x == 5 else 0 for x in responseMatrix[day]]
            nTooEarly[day] = sum(tooEarlies[day])
            missTrials[day] = [1 if x==0 else 0 for x in responseMatrix[day]]
            nMissTrial[day] = sum(missTrials[day])
        else: 
            hits[day] = np.nan; nHit[day] = np.nan; misses[day] = np.nan; nMiss[day] = np.nan; correctRejects[day] = np.nan; nCorrectReject[day] = np.nan; falseAlarms[day] = np.nan; nFalseAlarm[day] = np.nan; tooEarlies[day] = np.nan; nTooEarly[day] = np.nan; missTrials[day] = np.nan; nMissTrial[day] = np.nan;       
            
        tType[day] = allDataSheet.loc[day, 'trialType'] #1 go, 0 no/tooEarly (tooEarly if trialOutcome=5)
        if not 2 in tType[day]:
            if isinstance(responseMatrix[day], np.ndarray):
                for trial in range(nTrials[day]):
                    if tType[day][trial] == 1:
                        tType_calc[day] = np.append(tType_calc[day], 1)
                    elif tType[day][trial] == 0 and responseMatrix[day][trial] == 0:
                        tType_calc[day] = np.append(tType_calc[day], 0)
                    elif tType[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                        tType_calc[day] = np.append(tType_calc[day], 2)    
                    else:
                        tType_calc[day] = np.append(tType_calc[day], 0)
            else:
                tType_calc[day] = tType[day]
        elif 2 in tType[day]:
            tType_calc[day] = tType[day] #for ease
        tBeamBreak[day] = allDataSheet.loc[day, 'beamBreakTime'] #
        tTerminalBeamBreak[day] = allDataSheet.loc[day, 'terminalBeamBreakTime']
        
###### This is a few lines to separate tooEarlies that happen before the tone (between start cue and tone) and those that
    # happen immediately after (tone + 0.2). OOC right now, will need to be changed slightly when I have a variable that lines nose pokes
    # (tResponse) to cue,trial start, etc.
    """    
        stimuli = [1 if x>0 else 0 for x in rawTone]
        nStimuli = sum(stimuli)
        tooEarlies_afterStartCue = [
            1 if x == 1 and y == 1 else 0
            for x,y in zip(tooEarlies,stimuli)
            ]
        nTooEarly_afterStartCue = sum(tooEarlies_afterStartCue)
        tooEarlies_duringStartCue = [
            1 if x == 1 and y == 0 else 0
            for x,y in zip(tooEarlies,stimuli)
            ]
        nTooEarly_duringStartCue = sum(tooEarlies_duringStartCue)
    """
    
    #  time for some maths
    for day in range(len(allDataSheet)):    
        if 'ST' in str(allDataSheet.loc[day,'Experiment']):
            phaseIndex[day] = 1
            rateAccuracy[day] = np.append(rateAccuracy[day], nTrials[day]/maxTrials[day]) 
            rateAccuracyGo[day] = np.append(rateAccuracyGo[day], 0)
            rateAccuracyNo[day] = np.append(rateAccuracyNo[day], 0)
            rateTooEarly[day] = np.append(rateTooEarly[day], 0)
        elif 'GT' in str(allDataSheet.loc[day,'Experiment']):
            phaseIndex[day] = 2
            rateAccuracy[day] = np.append(rateAccuracy[day], np.count_nonzero(rawReward[day])/nTrials[day])
            rateAccuracyGo[day] = np.append(rateAccuracyGo[day], nHit[day]/(nHit[day]+nMiss[day]))
            rateAccuracyNo[day] = np.append(rateAccuracyNo[day], 0)
            rateTooEarly[day] = np.append(rateTooEarly[day], np.count_nonzero(tooEarlies[day])/nTrials[day])
        elif 'NT' in str(allDataSheet.loc[day,'Experiment']):
            phaseIndex[day] = 3
            rateAccuracy[day] = np.append(rateAccuracy[day], np.count_nonzero(rawReward[day])/nTrials[day]) 
            rateAccuracyGo[day] = np.append(rateAccuracyGo[day], nHit[day]/(nHit[day]+nMiss[day]))
            rateAccuracyNo[day] = np.append(rateAccuracyNo[day], nCorrectReject[day]/(nCorrectReject[day]+nFalseAlarm[day]))
            rateTooEarly[day] = np.append(rateTooEarly[day], np.count_nonzero(tooEarlies[day])/nTrials[day])
        elif 'RT' in str(allDataSheet.loc[day,'Experiment']):
            phaseIndex[day] = 4
            rateAccuracy[day] = np.append(rateAccuracy[day], np.count_nonzero(rawReward[day])/nTrials[day]) 
            rateAccuracyGo[day] = np.append(rateAccuracyGo[day], nHit[day]/(nHit[day]+nMiss[day]))
            rateAccuracyNo[day] = np.append(rateAccuracyNo[day], nCorrectReject[day]/(nCorrectReject[day]+nFalseAlarm[day]))
            rateTooEarly[day] = np.append(rateTooEarly[day], np.count_nonzero(tooEarlies[day])/nTrials[day])
        
######  time for some basic graphs (line for training trajectory)
    # Total accuracy (does not discriminate between go trials and no-go trials)
    plt.figure(figsize=(20, 12))
    plt.plot(list(range(1, len(allDataSheet)+1)), rateAccuracy, color='purple', linewidth=4, markersize=10, marker='o')
    plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
    
    plt.title('Accuracy across all trials training sessions'); plt.xlabel('Days'); plt.ylabel('% Correct')
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks(list(range(1, len(allDataSheet)+1))); ax.set_yticks(np.linspace(0,1,11)); ax.set_xlim(0, len(allDataSheet)); ax.set_ylim(0, 1)
    
    if 1 in phaseIndex:
        phaseST = plt.axvspan((1+phaseIndex.index(1)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(1))), color='yellow', alpha=0.1, label='ST')  #spout training
    if 2 in phaseIndex:
        phaseGT = plt.axvspan((1+phaseIndex.index(2)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2))), color='yellow', alpha=0.25, label='GT') #go training
    if 3 in phaseIndex:
        phaseNT = plt.axvspan((1+phaseIndex.index(3)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(3))), color='yellow', alpha=0.4, label='NT')  #no-go training
    if 4 in phaseIndex:
        phaseRT = plt.axvspan((1+phaseIndex.index(4)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))), color='yellow', alpha=0.55, label='RT') #reversal training
    
    plt.legend(fontsize=24,loc='upper left');
    
    plt.savefig('{directory}fullSessionAccuracy_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    # Go trial accuracy
    plt.figure(figsize=(20, 12))
    plt.plot(list(range(1, len(allDataSheet)+1)), rateAccuracyGo, color='blue', linewidth=4, markersize=10, marker='o', label='go')
    plt.plot(list(range(1, len(allDataSheet)+1)), rateAccuracyNo, color='red', linewidth=4, markersize=10, marker='o', label='no-go')
    plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
    
    plt.title('Accuracy across training sessions',fontsize=48); plt.xlabel('Days',fontsize=36); plt.ylabel('% Correct',fontsize=36)
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks(list(range(1, len(allDataSheet)+1))); ax.set_yticks(np.linspace(0,1,11)); ax.set_xlim(0, len(allDataSheet)); ax.set_ylim(0, 1)
    plt.tick_params(axis='both', labelsize=36);
    
    if 1 in phaseIndex:
        phaseST = plt.axvspan((1+phaseIndex.index(1)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(1))), color='yellow', alpha=0.1, label='ST')  #spout training
    if 2 in phaseIndex:
        phaseGT = plt.axvspan((1+phaseIndex.index(2)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2))), color='yellow', alpha=0.25, label='GT') #go training
    if 3 in phaseIndex:
        phaseNT = plt.axvspan((1+phaseIndex.index(3)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(3))), color='yellow', alpha=0.4, label='NT')  #no-go training
    if 4 in phaseIndex:
        phaseRT = plt.axvspan((1+phaseIndex.index(4)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))), color='yellow', alpha=0.55, label='RT') #reversal training
    
    plt.legend(fontsize=24,loc='upper left');
    
    plt.savefig('{directory}goNoGoTrialsAccuracy_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
          
    # Rate of too early responses  
    plt.figure(figsize=(20, 12))
    plt.plot(list(range(1, len(allDataSheet)+1)), rateTooEarly, color='black')
    
    plt.title('Percent too early across training sessions',fontsize=48); plt.xlabel('Days',fontsize=36); plt.ylabel('% Too Early',fontsize=36)
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks(list(range(1, len(allDataSheet)+1))); ax.set_yticks(np.linspace(0,1,11)); ax.set_xlim(0, len(allDataSheet)); ax.set_ylim(0, 1)
    plt.tick_params(axis='both', labelsize=36);
    
    if 1 in phaseIndex:
        phaseST = plt.axvspan((1+phaseIndex.index(1)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(1))), color='yellow', alpha=0.1, label='ST')  #spout training
    if 2 in phaseIndex:
        phaseGT = plt.axvspan((1+phaseIndex.index(2)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2))), color='yellow', alpha=0.25, label='GT') #go training
    if 3 in phaseIndex:
        phaseNT = plt.axvspan((1+phaseIndex.index(3)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(3))), color='yellow', alpha=0.4, label='NT')  #no-go training
    if 4 in phaseIndex:
        phaseRT = plt.axvspan((1+phaseIndex.index(4)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))), color='yellow', alpha=0.55, label='RT') #reversal training
    
    plt.legend(fontsize=24,loc='upper left');
    
    plt.savefig('{directory}fullSessionTooEarly_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
###### Now some scatter plots to align events to trial start
    
    tEnd_calc = [[] for _ in range(len(allDataSheet))];
    rawBeamBreak_calc = [[] for _ in range(len(allDataSheet))]; tBeamBreak_calc = [[] for _ in range(len(allDataSheet))]; tBeamBreak_calcIndex = [[] for _ in range(len(allDataSheet))];
    
    for day in range(len(allDataSheet)):
        if any(x<0 for x in tReward[day]):
            for trial in range(nTrials[day]):
                if tReward[day][trial] == 0 or rawEnd[day][trial] == 3600:
                    tEnd_calc[day] = np.append(tEnd_calc[day], -999) #need to format unwanted data so that it preserves trial structure (nan does not work for scatter plots)
                elif tReward[day][trial] > 0: #reward delivered
                    tEnd_calc[day] = np.append(tEnd_calc[day], tReward[day][trial]+0.01) #if reward was released, tEnd is reward + 500 msec (in line with code update on 240102 which added a beam break after 0.5s of reward to initiate trial end)
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 1: #no reward during go trial
                    tEnd_calc[day] = np.append(tEnd_calc[day], tStart[day][trial]+7.5) #if no reward, and go trial, trial ends 7.5 seconds [start > 0.5 > trialLight on > 2 > tone on > window close] after start
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 2: #no reward during nogo trial
                    tEnd_calc[day] = np.append(tEnd_calc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #it is written with 1 sec after failed trial to allow remaining formulas to finish (due to overlapping logicals that do not interact)
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 0 and responseMatrix[day][trial] == 5: #too early         ^
                    tEnd_calc[day] = np.append(tEnd_calc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #same as above |
        
    for day in range(len(allDataSheet)):
        if isinstance(tBeamBreak[day], np.ndarray):
            if 'ST' in str(allDataSheet.loc[day,'Experiment']):
                for trial in range(nTrials[day]):
                    afterStart = tBeamBreak[day]>rawStart[day][trial] #reward delivered at start
                    firstResp = next((resp for resp, log in enumerate(afterStart) if log), None)
                    if firstResp == None:
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)
                    else:
                        if not rawStart[day][trial] == 0:
                            rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day], tBeamBreak[day][firstResp])
                            tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day], rawBeamBreak_calc[day][trial]-rawStart[day][trial])
                        else:
                            rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                            tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999) 
            else:
                for trial in range(nTrials[day]):
                    afterStart = tBeamBreak[day]>rawStart[day][trial]
                    beforeEnd = tBeamBreak[day]<(rawStart[day][trial]+tEnd_calc[day][trial])
                    if not np.logical_and(afterStart, beforeEnd).any():
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)
                    else:
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day], tBeamBreak[day][np.logical_and(afterStart, beforeEnd)][0])
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day], rawBeamBreak_calc[day][trial]-rawStart[day][trial])           
        elif np.isnan(tBeamBreak[day]):
            for trial in range(nTrials[day]):
                tBeamBreak_calc[day]= np.append(tBeamBreak_calc[day], -999)
    
    """
    tTone_calc = [[] for _ in range(len(allDataSheet))]; tReward_calc = [[] for _ in range(len(allDataSheet))]; tTerminalBeamBreak_calc = [[] for _ in range(len(allDataSheet))]; 
    tEnd_goCalc = [[] for _ in range(len(allDataSheet))]; tTone_goCalc = [[] for _ in range(len(allDataSheet))]; tReward_goCalc = [[] for _ in range(len(allDataSheet))]; tTerminalBeamBreak_goCalc = [[] for _ in range(len(allDataSheet))]; 
    tEnd_noCalc = [[] for _ in range(len(allDataSheet))]; tTone_noCalc = [[] for _ in range(len(allDataSheet))]; tReward_noCalc = [[] for _ in range(len(allDataSheet))]; tTerminalBeamBreak_noCalc = [[] for _ in range(len(allDataSheet))]; 
    
    # All trial events
    for day in range(len(allDataSheet)):
        # start of trial
        startIcon = plt.scatter(tStart[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawStart[day])+1)), marker='|', color='black', alpha=0.5, label='trial start')
        # end of trial
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tReward[day][trial] == 0:
                    tEnd_calc[day] = np.append(tEnd_calc[day], -999) #need to format unwanted data so that it preserves trial structure (nan does not work for scatter plots)
                elif tReward[day][trial] > 0: #reward delivered
                    tEnd_calc[day] = np.append(tEnd_calc[day], tReward[day][trial]+0.01) #if reward was released, tEnd is reward + 500 msec (in line with code update on 240102 which added a beam break after 0.5s of reward to initiate trial end)
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 1: #no reward during go trial
                    tEnd_calc[day] = np.append(tEnd_calc[day], tStart[day][trial]+7.5) #if no reward, and go trial, trial ends 7.5 seconds [start > 0.5 > trialLight on > 2 > tone on > window close] after start
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 2: #no reward during nogo trial
                    tEnd_calc[day] = np.append(tEnd_calc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #it is written with 1 sec after failed trial to allow remaining formulas to finish (due to overlapping logicals that do not interact)
                elif tReward[day][trial] < 0 and tType_calc[day][trial] == 0 and responseMatrix[day][trial] == 5: #too early         ^
                    tEnd_calc[day] = np.append(tEnd_calc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #same as above |
            endIcon = plt.scatter(tEnd_calc[day][0:np.count_nonzero(tEnd_calc[day])], list(range(1, np.count_nonzero(tEnd_calc[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
        else:
            endIcon = plt.scatter(tEnd[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(tEnd[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
            
        # tone plays
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tTone[day][trial] == 0:
                    tTone_calc[day] = np.append(tTone_calc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                elif tTone[day][trial] < 0:
                    tTone_calc[day] = np.append(tTone_calc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                else:    
                    tTone_calc[day] = np.append(tTone_calc[day], tTone[day][trial])
        else:
            tTone_calc[day] = np.append(tTone_calc[day], tTone[day])
        toneIcon = plt.scatter(tTone_calc[day][0:np.count_nonzero(tTone_calc[day])], list(range(1, np.count_nonzero(tTone_calc[day])+1)), marker='d', color='purple', s=15, label='tone') # black diamond
    
        # beam break
        if isinstance(tBeamBreak[day], np.ndarray):
            for trial in range(maxTrials[day]):
                afterStart = tBeamBreak[day]>rawStart[day][trial]
                beforeEnd = tBeamBreak[day]<(rawStart[day][trial]+tEnd_calc[day][trial])
                if not np.logical_and(afterStart, beforeEnd).any():
                    rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                    tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)
                else:
                    rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day], tBeamBreak[day][np.logical_and(afterStart, beforeEnd)][0])
                    tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day], rawBeamBreak_calc[day][trial]-rawStart[day][trial])           
        elif np.isnan(tBeamBreak[day]):
            for trial in range(maxTrials[day]):
                tBeamBreak_calc[day]= np.append(tBeamBreak_calc[day], -999)
        beamIcon = plt.scatter(tBeamBreak_calc[day], list(range(1, np.count_nonzero(tBeamBreak_calc[day])+1)), marker='.', color='black', s=15, label='response') # black dot
            
        # reward
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tReward[day][trial] == 0:
                    tReward_calc[day] = np.append(tReward_calc[day], -999) 
                elif tReward[day][trial] < 0:
                    tReward_calc[day] = np.append(tReward_calc[day], -999) 
                else:    
                    tReward_calc[day] = np.append(tReward_calc[day], tReward[day][trial])
            rewardIcon = plt.scatter(tReward_calc[day], list(range(1, np.count_nonzero(tReward_calc[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle
        else:
            rewardIcon = plt.scatter(tReward[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawStart[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle        
            
        #plt.legend(loc='lower right')
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_ylim(0, (np.count_nonzero(rawStart[day])+1)); ax.set_xlim(-1, 9); ax.set_xticks([-1,0,1,2,3,4,5,6,7,8,9]); ax.set_ylabel('Trials'); ax.set_xlabel('Time from trial start (sec)')
        plt.title('rasterArray_allTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
        
        figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
        figName = '{directory}rasterArray_allTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
        if os.path.isfile(os.path.join(outputPath,figName)):
            continue
        else:    
            plt.savefig('{directory}rasterArray_allTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
        plt.show(); plt.close()
        
    # Go trial events
    for day in range(len(allDataSheet)):
        # start of trial
        startIcon = plt.scatter(tStart[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawStart[day])+1)), marker='|', color='black', alpha=0.5, label='trial start')
        # end of trial
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 1:
                    if tReward[day][trial] == 0:
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999)
                    elif tReward[day][trial] > 0 and tType_calc[day][trial] == 1: #reward delivered during go trial
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], tReward[day][trial]+0.01) #if reward was released, tEnd is reward + 500 msec (in line with code update on 240102 which added a beam break after 0.5s of reward to initiate trial end)
                    elif tReward[day][trial] > 0 and tType_calc[day][trial] == 2: #reward during no go trial
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999)
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 1: #no reward during go trial
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], tStart[day][trial]+7.5) #if no reward, and go trial, trial ends 7.5 seconds [start > 0.5 > trialLight on > 2 > tone on > window close] after start
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 2: #no reward during nogo trial
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999) #nogo trials don't count here
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 0 and responseMatrix[day][trial] == 5: #too early         
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999) #not bothering with too earlies when looking at only go trials
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999)
            endIcon = plt.scatter(tEnd_goCalc[day][0:np.count_nonzero(tEnd_goCalc[day])], list(range(1, np.count_nonzero(tEnd_goCalc[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
        else:
            endIcon = plt.scatter(tEnd[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(tEnd[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
            
        # tone plays
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 1:
                    if tTone[day][trial] == 0:
                        tTone_goCalc[day] = np.append(tTone_goCalc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                    elif tTone[day][trial] < 0:
                        tTone_goCalc[day] = np.append(tTone_goCalc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                    else:    
                        tTone_goCalc[day] = np.append(tTone_goCalc[day], tTone[day][trial])
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tTone_goCalc[day] = np.append(tTone_goCalc[day], -999)
        else:
            tTone_goCalc[day] = np.append(tTone_goCalc[day], tTone[day])            
        toneIcon = plt.scatter(tTone_goCalc[day][0:np.count_nonzero(tTone_calc[day])], list(range(1, np.count_nonzero(tTone_goCalc[day])+1)), marker='d', color='blue', s=15, label='tone') # black diamond
    
        # beam break
        if isinstance(tBeamBreak[day], np.ndarray):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 1:
                    afterStart = tBeamBreak[day]>rawStart[day][trial]
                    beforeEnd = tBeamBreak[day]<(rawStart[day][trial]+tEnd_calc[day][trial])
                    if not np.logical_and(afterStart, beforeEnd).any():
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)
                    else:
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day], tBeamBreak[day][np.logical_and(afterStart, beforeEnd)][0])
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day], rawBeamBreak_calc[day][trial]-rawStart[day][trial]) 
                elif tType_calc[day][trial] == 2 or tType_calc[day][trial] == 0:
                    tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)  
        elif np.isnan(tBeamBreak[day]):
            for trial in range(maxTrials[day]):
                tBeamBreak_calc[day]= np.append(tBeamBreak_calc[day], -999)
        beamIcon = plt.scatter(tBeamBreak_calc[day], list(range(1, np.count_nonzero(tBeamBreak_calc[day])+1)), marker='.', color='black', s=15, label='response') # black dot
    
        # reward
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 1:
                    if tReward[day][trial] == 0:
                        tReward_goCalc[day] = np.append(tReward_goCalc[day], -999) 
                    elif tReward[day][trial] < 0:
                        tReward_goCalc[day] = np.append(tReward_goCalc[day], -999) 
                    else:    
                        tReward_goCalc[day] = np.append(tReward_goCalc[day], tReward[day][trial])
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tReward_goCalc[day] = np.append(tReward_goCalc[day], -999)
            rewardIcon = plt.scatter(tReward_goCalc[day], list(range(1, np.count_nonzero(tReward_goCalc[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle
        else:
            rewardIcon = plt.scatter(tReward[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawReward[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle        
            #for trial in range(maxTrials[day]):
            #    tReward_goCalc[day] = np.append(tReward_goCalc[day], -999)
            #rewardIcon = plt.scatter(tReward_goCalc[day], list(range(1, np.count_nonzero(tReward_goCalc[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle
            
        #plt.legend(loc='lower right')
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_ylim(0, len(tEnd_goCalc[day])+1); ax.set_xlim(-1, 9); ax.set_xticks([-1,0,1,2,3,4,5,6,7,8,9]); ax.set_ylabel('Trials'); ax.set_xlabel('Time from trial start (sec)')
        plt.title('rasterArray_goTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
        
        figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
        figName = '{directory}rasterArray_goTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
        if os.path.isfile(os.path.join(outputPath,figName)):
            continue
        else:    
            plt.savefig('{directory}rasterArray_goTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
        plt.show(); plt.close()
    
    # No-go trial events
    for day in range(len(allDataSheet)):
        # start of trial
        startIcon = plt.scatter(tStart[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawStart[day])+1)), marker='|', color='black', alpha=0.5, label='trial start')
        # end of trial
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    if tReward[day][trial] == 0:
                        tEnd_noCalc[day] = np.append(tEnd_noCalc[day], -999)
                    elif tReward[day][trial] > 0 and tType_calc[day][trial] == 2: #reward delivered for no go trial
                        tEnd_noCalc[day] = np.append(tEnd_noCalc[day], tReward[day][trial]+0.01) #if reward was released, tEnd is reward + 500 msec (in line with code update on 240102 which added a beam break after 0.5s of reward to initiate trial end)
                    elif tReward[day][trial] > 0 and tType_calc[day][trial] == 1: #reward during no go trial
                        tEnd_goCalc[day] = np.append(tEnd_goCalc[day], -999)
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 1: #no reward during go trial
                        tEnd_noCalc[day] = np.append(tEnd_noCalc[day], -999) #go trials don't count here
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 2: #no reward during nogo trial
                        tEnd_noCalc[day] = np.append(tEnd_noCalc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #it is written with 1 sec after failed trial to allow remaining formulas to finish (due to overlapping logicals that do not interact)
                    elif tReward[day][trial] < 0 and tType_calc[day][trial] == 0 and responseMatrix[day][trial] == 5: #too early             ^
                        tEnd_noCalc[day] = np.append(tEnd_noCalc[day], tTerminalBeamBreak[day][trial]-rawStart[day][trial]+1) #same as above |
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tTone_noCalc[day] = np.append(tTone_noCalc[day], -999)            
            endIcon = plt.scatter(tEnd_noCalc[day][0:np.count_nonzero(tEnd_noCalc[day])], list(range(1, np.count_nonzero(tEnd_noCalc[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
        else:
            endIcon = plt.scatter(tEnd[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(tEnd[day])+1)), marker='|', color='black', alpha=0.5, label='trial end') # black |
            
        # tone plays
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    if tTone[day][trial] == 0:
                        tTone_noCalc[day] = np.append(tTone_noCalc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                    elif tTone[day][trial] < 0:
                        tTone_noCalc[day] = np.append(tTone_noCalc[day], -999) #negative values (i.e. no reward delivered) will fall out of range after scaling the y-axis to -1:9 sec
                    else:    
                        tTone_noCalc[day] = np.append(tTone_noCalc[day], tTone[day][trial])
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tTone_noCalc[day] = np.append(tTone_noCalc[day], -999)        
        else:
            tTone_noCalc[day] = np.append(tTone_noCalc[day], tTone[day]) 
        toneIcon = plt.scatter(tTone_noCalc[day][0:np.count_nonzero(tTone_calc[day])], list(range(1, np.count_nonzero(tTone_noCalc[day])+1)), marker='d', color='red', s=15, label='tone') # black diamond
    
        # beam break
        if isinstance(tBeamBreak[day], np.ndarray):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 2:
                    afterStart = tBeamBreak[day]>rawStart[day][trial]
                    beforeEnd = tBeamBreak[day]<(rawStart[day][trial]+tEnd_calc[day][trial])
                    if not np.logical_and(afterStart, beforeEnd).any():
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day],-999) #no response
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)
                    else:
                        rawBeamBreak_calc[day] = np.append(rawBeamBreak_calc[day], tBeamBreak[day][np.logical_and(afterStart, beforeEnd)][0])
                        tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day], rawBeamBreak_calc[day][trial]-rawStart[day][trial]) 
                elif tType_calc[day][trial] == 1 or tType_calc[day][trial] == 0:
                    tBeamBreak_calc[day] = np.append(tBeamBreak_calc[day],-999)  
        elif np.isnan(tBeamBreak[day]):
            for trial in range(maxTrials[day]):
                tBeamBreak_calc[day]= np.append(tBeamBreak_calc[day], -999)
        beamIcon = plt.scatter(tBeamBreak_calc[day], list(range(1, np.count_nonzero(tBeamBreak_calc[day])+1)), marker='.', color='black', s=15, label='response') # black dot
    
        # reward
        if any(x<0 for x in tReward[day]):
            for trial in range(maxTrials[day]):
                if tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    if tReward[day][trial] == 0:
                        tReward_noCalc[day] = np.append(tReward_noCalc[day], -999) 
                    elif tReward[day][trial] < 0:
                        tReward_noCalc[day] = np.append(tReward_noCalc[day], -999) 
                    else:    
                        tReward_noCalc[day] = np.append(tReward_noCalc[day], tReward[day][trial])
                elif tType_calc[day][trial] == 0 and not responseMatrix[day][trial] == 5:
                    tTone_noCalc[day] = np.append(tTone_noCalc[day], -999)        
            rewardIcon = plt.scatter(tReward_noCalc[day], list(range(1, np.count_nonzero(tReward_noCalc[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle
        else:
            rewardIcon = plt.scatter(tReward[day][0:np.count_nonzero(rawStart[day])], list(range(1, np.count_nonzero(rawReward[day])+1)), marker='1', color='yellow', s=15, label='reward') # blue triangle        
            
        #plt.legend(loc='lower right')
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_ylim(0, len(tEnd_noCalc[day])+1); ax.set_xlim(-1, 9); ax.set_xticks([-1,0,1,2,3,4,5,6,7,8,9]); ax.set_ylabel('Trials'); ax.set_xlabel('Time from trial start (sec)')
        plt.title('rasterArray_nogoTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
        
        figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
        figName = '{directory}rasterArray_nogoTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
        if os.path.isfile(os.path.join(outputPath,figName)):
            continue
        else:    
            plt.savefig('{directory}rasterArray_nogoTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
        plt.show(); plt.close()    
    """
    
###### Here I'm going to plot a PSTH of these response frequencies, so nTrials (not including tooEarly) bin into 500msec bits and averaged
    
    psthBin = 200 #msec
    scaleFactor = 100 #binning in msec scale
    psthBinArray = np.arange(-10,91,psthBin/scaleFactor) # 
    psthBinStart = psthBinArray[0]
    psthBinEnd = psthBinArray[::-1][0]
    
    psthArray_resp = [[] for _ in range(len(allDataSheet))]; 
    psthArray_index = [[] for _ in range(len(allDataSheet))]; 
    psthArray_beamBreakParse = [[] for _ in range(len(allDataSheet))]; 
    psthArray_trialParse = [[] for _ in range(len(allDataSheet))]; 
    for day in range(len(allDataSheet)):
        psthArray_resp[day] = [[] for _ in range(nTrials[day])]  
        psthArray_index[day] = [[] for _ in range(nTrials[day])]  
        psthArray_beamBreakParse[day] = [[] for _ in range(nTrials[day])]  
        psthArray_trialParse[day] = [[] for _ in range(nTrials[day])]  
        for trial in range(nTrials[day]):
            if tBeamBreak_calc[day][trial] > 0:
                psthArray_beamBreakParse[day] = np.append(psthArray_beamBreakParse[day], tBeamBreak_calc[day][trial])
                psthArray_trialParse[day] = np.append(psthArray_trialParse[day], trial)
                psthArray_resp[day] = np.append(psthArray_resp[day], math.floor(tBeamBreak_calc[day][trial]*10))
                if tBeamBreak_calc[day][trial] < 0.5:
                    psthArray_index[day] = np.append(psthArray_index[day], 0)
                else:
                    psthArray_index[day] = np.append(psthArray_index[day], tType_calc[day][trial])
            elif tBeamBreak_calc[day][trial] == -999:
                continue
            
    psthArray_resp2D = [np.zeros((len(psthArray_resp[day]),len(psthBinArray))) for day in range(len(allDataSheet))]  
    psthArray_respBinned2D = [[] for _ in range(len(allDataSheet))]; psthArray_respBinned2D_go = [[] for _ in range(len(allDataSheet))]; psthArray_respBinned2D_no = [[] for _ in range(len(allDataSheet))]  
    for day in range(len(allDataSheet)):    
        if isinstance(psthArray_resp[day], np.ndarray):
            for trial in range(len(psthArray_resp[day])):
                for bins, binStart in enumerate(psthBinArray):
                    if binStart <= psthArray_resp[day][trial] <= binStart + int(psthBin/scaleFactor):
                        psthArray_resp2D[day][trial,bins] = 1 
                        break
        psthArray_respBinned2D[day] = np.sum(psthArray_resp2D[day], axis=0)
        psthArray_respBinned2D_go[day] = np.sum(psthArray_resp2D[day][psthArray_index[day]==1], axis=0)
        psthArray_respBinned2D_no[day] = np.sum(psthArray_resp2D[day][psthArray_index[day]==2], axis=0)
        
    psthCount_resp = [[] for _ in range(len(allDataSheet))]; psthAverage_resp = [[] for _ in range(len(allDataSheet))]
    for day in range(len(allDataSheet)):  
        for bin in psthBinArray:
            psthCount_resp[day] = np.append(psthCount_resp[day], np.count_nonzero(psthArray_resp[day]==bin))
            #psthAverage_resp[day] = np.append(psthAverage_resp[day], psthCount_resp[day][bin]/nTrials[day])   
    
###### time to plot the response histogram
    for day in range(len(allDataSheet)):
        
        if isinstance(psthArray_resp[day], np.ndarray):
            figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
            figName = '{directory}psthFirstResponse_allTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
            if os.path.isfile(os.path.join(outputPath,figName)):
                continue
            else:
                #all trials                
                plt.figure(figsize=(10, 6))
                plt.bar(np.linspace(-1,9,51), psthArray_respBinned2D[day], width=(psthBin/scaleFactor)/10, color='purple', edgecolor='purple', align='edge')
                plt.axvline(0, color='black', linestyle='-', linewidth=2)
                plt.axvline(2.5, color='black', linestyle='--', linewidth=2)
            
                ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
                ax.set_ylim(0, np.max(psthArray_respBinned2D[day])+1); ax.set_yticks(np.arange(0,np.max(psthArray_respBinned2D[day])+1)); ax.set_ylabel('# responses');
                ax.set_xlim(-1.0,9.0); ax.set_xticks(np.arange(psthBinStart/10,(psthBinEnd/10)+1)); ax.set_xlabel('Time from trial start (sec)');
                plt.title('psthFirstResponse_allTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
                    
                plt.savefig('{directory}psthFirstResponse_allTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
                plt.show() 
                plt.close()
        
            figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
            figName = '{directory}psthFirstResponse_goTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
            if os.path.isfile(os.path.join(outputPath,figName)):
                continue
            else:
                # go only
                plt.figure(figsize=(10, 6))
                plt.bar(np.linspace(-1,9,51), psthArray_respBinned2D_go[day], width=(psthBin/scaleFactor)/10, color='blue', edgecolor='blue', align='edge')
                plt.axvline(0, color='black', linestyle='-', linewidth=2)
                plt.axvline(2.5, color='black', linestyle='--', linewidth=2)
            
                ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
                ax.set_ylim(0, np.max(psthArray_respBinned2D[day])+1); ax.set_yticks(np.arange(0,np.max(psthArray_respBinned2D[day])+1)); ax.set_ylabel('# responses');
                ax.set_xlim(-1.0,9.0); ax.set_xticks(np.arange(psthBinStart/10,(psthBinEnd/10)+1)); ax.set_xlabel('Time from trial start (sec)');
                plt.title('psthFirstResponse_goTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
           
                plt.savefig('{directory}psthFirstResponse_goTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)    
                plt.show()
                plt.close()
        
            figFiles = [file for file in os.listdir(outputPath) if file.endswith(tuple('.pdf'))]
            figName = '{directory}psthFirstResponse_noTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse)
            if os.path.isfile(os.path.join(outputPath,figName)):
                continue
            else: 
                # no only
                plt.figure(figsize=(10, 6))
                plt.bar(np.linspace(-1,9,51), psthArray_respBinned2D_no[day], width=(psthBin/scaleFactor)/10, color='red', edgecolor='red', align='edge')
                plt.axvline(0, color='black', linestyle='-', linewidth=2)
                plt.axvline(2.5, color='black', linestyle='--', linewidth=2)
            
                ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
                ax.set_ylim(0, np.max(psthArray_respBinned2D[day])+1); ax.set_yticks(np.arange(0,np.max(psthArray_respBinned2D[day])+1)); ax.set_ylabel('# responses');
                ax.set_xlim(-1.0,9.0); ax.set_xticks(np.arange(psthBinStart/10,(psthBinEnd/10)+1)); ax.set_xlabel('Time from trial start (sec)');
                plt.title('psthFirstResponse_noTrials_{date}'.format(date=str(allDataSheet.loc[day,'Date'])))
           
                plt.savefig('{directory}psthFirstResponse_noTrials_{date}_{experiment}{mouse}.pdf'.format(directory=outputPath, date=str(allDataSheet.loc[day,'Date']), experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)   
                plt.show()
                plt.close()
    
###### time to plot a line graph to show the change in response number across these sessions (similar to accuracy but different)
    psthArray_summedResp_go = [[] for _ in range(len(allDataSheet))]; psthArray_summedResp_no = [[] for _ in range(len(allDataSheet))]; psthArray_summedResp = [[] for _ in range(len(allDataSheet))]
    for day in range(len(allDataSheet)):
        psthArray_summedResp[day] = np.append(psthArray_summedResp[day], np.sum(psthArray_respBinned2D[day]))
        psthArray_summedResp_go[day] = np.append(psthArray_summedResp_go[day], np.sum(psthArray_respBinned2D_go[day]))
        psthArray_summedResp_no[day] = np.append(psthArray_summedResp_no[day], np.sum(psthArray_respBinned2D_no[day]))
        
        #plot
    plt.figure(figsize=(10, 6))
    plt.plot(list(range(1, len(allDataSheet)+1)), psthArray_summedResp_go, color='blue', label='go')
    plt.plot(list(range(1, len(allDataSheet)+1)), psthArray_summedResp_no, color='red', label='no')
    
    plt.title('Number of responses across training sessions'); plt.xlabel('Days'); plt.ylabel('# Responses')
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xlim(0, len(allDataSheet)); ax.set_ylim(0, 180);
    ax.set_xticks(list(range(1, len(allDataSheet)+1))); 
    
    if 1 in phaseIndex:
        phaseST = plt.axvspan((1+phaseIndex.index(1)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(1))), color='yellow', alpha=0.1, label='ST')  #spout training
    if 2 in phaseIndex:
        phaseGT = plt.axvspan((1+phaseIndex.index(2)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2))), color='yellow', alpha=0.25, label='GT') #go training
    if 3 in phaseIndex:
        phaseNT = plt.axvspan((1+phaseIndex.index(3)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(3))), color='yellow', alpha=0.4, label='NT')  #no-go training
    if 4 in phaseIndex:
        phaseRT = plt.axvspan((1+phaseIndex.index(4)), (1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))), color='yellow', alpha=0.55, label='RT') #reversal training
    
    plt.legend(loc='upper left');
    
    plt.savefig('{directory}fullSessionResponseCount_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close() 
    
    sSegments = 30;   
    iSegment = [[[sS + iT - 1 for iT in range(sSegments) if sS+iT < nTrials[d]] for sS in range(1, (nTrials[d]), sSegments)] for d in range(len(allDataSheet))]
    rateAccuracy_segment = [[[np.count_nonzero(rawReward[d][iSegment[d][s]])/len(iSegment[d][s])] for s in range(len(iSegment[d]))] for d in range(len(allDataSheet))]
    rateAccuracy_segment_percentChange = [[[np.round(rateAccuracy_segment[d][s][0]-rateAccuracy_segment[d][0][0],3)] for s in range(len(rateAccuracy_segment[d]))] for d in range(len(allDataSheet))]
    
    uniqVal, iFirst = np.unique(phaseIndex, return_index=True)
                                    
   
    plt.figure(figsize=(7, 6))
    for i in range(len(np.linspace(1,10,10))):
        plt.axvline(i, color='0.85', linestyle='-', linewidth=2)
    for d in range(len(allDataSheet)):
        if phaseIndex[d] > 1: # skip ST with phaseIndex>1
            if d in iFirst:
                plt.plot(rateAccuracy_segment[d], label=d)
        else:
            continue
    plt.legend(title='Session', loc='upper right', bbox_to_anchor=(1.1, 1));
    plt.xlabel('block (20 trials)')
    plt.ylabel('accuracy (%)')
    plt.title('From first session of each phase : {mouse}'.format(mouse=thisMouse))   
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xlim(0, 6); ax.set_ylim(0, 1);
    ax.set_yticks(np.linspace(0, 1, 11)); 
    plt.savefig('{directory}firstSessionAccuracyBinned_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
   
    #accuracy segments across go days
    
    if 2 in phaseIndex:
        nGo = np.round(np.linspace(phaseIndex.index(2),1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2)),(1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(2))-phaseIndex.index(2))+1)).astype(int)
        plt.figure(figsize=(7, 6))
        plt.xlabel('block (20 trials)')
        plt.ylabel('accuracy (%)')
        plt.title('From go phase : {mouse}'.format(mouse=thisMouse))   
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xlim(0, 6); ax.set_ylim(0, 1);
        ax.set_yticks(np.linspace(0, 1, 11)); 
        for i in range(len(np.linspace(1,10,10))):
            plt.axvline(i, color='0.85', linestyle='-', linewidth=2)
        for d in range(len(nGo)):
            plt.plot(rateAccuracy_segment[nGo[d]], label=nGo[d])
    
        plt.legend(title='Session', loc='upper right', bbox_to_anchor=(1.1, 1));
        plt.savefig('{directory}goSessionAccuracyBinned_{experiment}{mouse}.pdf'.format(directory=outputPath, experiment=str(allDataSheet.loc[0,'Experiment'])[0:str(allDataSheet.loc[0,'Experiment']).index('-')], mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        
    #accuracy segments from first 6 days of no-go phase
    
    if 3 in phaseIndex:
        nNoGo = np.round(np.linspace(phaseIndex.index(3),(len(phaseIndex) - 1 - phaseIndex[::-1].index(3)),(1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(3))-phaseIndex.index(3)))).astype(int)
        cMap=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']       
        
        rateAccuracyGo_segment = [[] for _ in range(len(nNoGo))]; rateAccuracyNo_segment = [[] for _ in range(len(nNoGo))]
        for d in range(len(nNoGo)):
                for s in range(len(iSegment[nNoGo[d]])):
                    if (sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==1)+sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==3))>0:
                        rateAccuracyGo_segment[d].append(sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==1)/(sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==1)+sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==3)))
                    if (sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==2)+sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==4))>0:    
                        rateAccuracyNo_segment[d].append(sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==2)/(sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==2)+sum(responseMatrix[nNoGo[d]][iSegment[nNoGo[d]][s]]==4)))
                                           
        fig,axs=plt.subplots(2,3,figsize=(10, 10))
        
        #
        if len(nNoGo)>=1:
            axs[0,0].plot(rateAccuracyGo_segment[0], color=cMap[0], label=nNoGo[0]+1)
            axs[0,0].plot(rateAccuracyNo_segment[0], color=cMap[0], linestyle='dashed')
            axs[0,0].set_ylabel('accuracy (%)')
            axs[0,0].spines['top'].set_visible(False); axs[0,0].spines['right'].set_visible(False);
            axs[0,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,0].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nNoGo)>=2:
            axs[0,1].plot(rateAccuracyGo_segment[1], color=cMap[1], label=nNoGo[1]+1)
            axs[0,1].plot(rateAccuracyNo_segment[1], color=cMap[1], linestyle='dashed')
            axs[0,1].spines['top'].set_visible(False); axs[0,1].spines['right'].set_visible(False);
            axs[0,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,1].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nNoGo)>=3:
            axs[0,2].plot(rateAccuracyGo_segment[2], color=cMap[2], label=nNoGo[2]+1)
            axs[0,2].plot(rateAccuracyNo_segment[2], color=cMap[2], linestyle='dashed')
            axs[0,2].spines['top'].set_visible(False); axs[0,2].spines['right'].set_visible(False);
            axs[0,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,2].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nNoGo)>=4:
            axs[1,0].plot(rateAccuracyGo_segment[3], color=cMap[3], label=nNoGo[3]+1)
            axs[1,0].plot(rateAccuracyNo_segment[3], color=cMap[3], linestyle='dashed')
            axs[1,0].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,0].set_ylabel('accuracy (%)')
            axs[1,0].spines['top'].set_visible(False); axs[1,0].spines['right'].set_visible(False);
            axs[1,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,0].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nNoGo)>=5:
            axs[1,1].plot(rateAccuracyGo_segment[4], color=cMap[4], label=nNoGo[4]+1)
            axs[1,1].plot(rateAccuracyNo_segment[4], color=cMap[4], linestyle='dashed')
            axs[1,1].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,1].spines['top'].set_visible(False); axs[1,1].spines['right'].set_visible(False);
            axs[1,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,1].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nNoGo)>=6:
            axs[1,2].plot(rateAccuracyGo_segment[5], color=cMap[5], label=nNoGo[5]+1)
            axs[1,2].plot(rateAccuracyNo_segment[5], color=cMap[5], linestyle='dashed')
            axs[1,2].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,2].spines['top'].set_visible(False); axs[1,2].spines['right'].set_visible(False);
            axs[1,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,2].set_yticks(np.linspace(0, 1, 11)); 
                     
        fig.suptitle('From no-go phase : {mouse}'.format(mouse=thisMouse))   
        plt.savefig('{directory}noGoPhaseAccuracyBinned_firstSix_{mouse}.pdf'.format(directory=outputPath, mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
    
    
    #accuracy segments from first 6 days of reversal phase
    if 4 in phaseIndex:
        nRev = np.round(np.linspace(phaseIndex.index(4),(len(phaseIndex) - 1 - phaseIndex[::-1].index(4)),(1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))-phaseIndex.index(4)))).astype(int)
        cMap=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']       
        
        rateAccuracyGo_segment = [[] for _ in range(len(nRev))]; rateAccuracyNo_segment = [[] for _ in range(len(nRev))]
        for d in range(len(nRev)):
                for s in range(len(iSegment[nRev[d]])):
                    if sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==3)>0:
                        rateAccuracyGo_segment[d].append(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)/(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==3)))
                    if sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==4)>0:
                        rateAccuracyNo_segment[d].append(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)/(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==4)))
                                           
        fig,axs=plt.subplots(2,3,figsize=(10, 10))
        
        #
        if len(nRev)>=1:
            axs[0,0].plot(rateAccuracyGo_segment[0], color=cMap[0], label=nRev[0]+1)
            axs[0,0].plot(rateAccuracyNo_segment[0], color=cMap[0], linestyle='dashed')
            axs[0,0].set_ylabel('accuracy (%)')
            axs[0,0].spines['top'].set_visible(False); axs[0,0].spines['right'].set_visible(False);
            axs[0,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,0].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nRev)>=2:
            axs[0,1].plot(rateAccuracyGo_segment[1], color=cMap[1], label=nRev[1]+1)
            axs[0,1].plot(rateAccuracyNo_segment[1], color=cMap[1], linestyle='dashed')
            axs[0,1].spines['top'].set_visible(False); axs[0,1].spines['right'].set_visible(False);
            axs[0,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,1].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nRev)>=3:
            axs[0,2].plot(rateAccuracyGo_segment[2], color=cMap[2], label=nRev[2]+1)
            axs[0,2].plot(rateAccuracyNo_segment[2], color=cMap[2], linestyle='dashed')
            axs[0,2].spines['top'].set_visible(False); axs[0,2].spines['right'].set_visible(False);
            axs[0,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[0,2].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nRev)>=4:
            axs[1,0].plot(rateAccuracyGo_segment[3], color=cMap[3], label=nRev[3]+1)
            axs[1,0].plot(rateAccuracyNo_segment[3], color=cMap[3], linestyle='dashed')
            axs[1,0].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,0].set_ylabel('accuracy (%)')
            axs[1,0].spines['top'].set_visible(False); axs[1,0].spines['right'].set_visible(False);
            axs[1,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,0].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nRev)>=5:
            axs[1,1].plot(rateAccuracyGo_segment[4], color=cMap[4], label=nRev[4]+1)
            axs[1,1].plot(rateAccuracyNo_segment[4], color=cMap[4], linestyle='dashed')
            axs[1,1].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,1].spines['top'].set_visible(False); axs[1,1].spines['right'].set_visible(False);
            axs[1,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,1].set_yticks(np.linspace(0, 1, 11)); 
    
        #
        if len(nRev)>=6:
            axs[1,2].plot(rateAccuracyGo_segment[5], color=cMap[5], label=nRev[5]+1)
            axs[1,2].plot(rateAccuracyNo_segment[5], color=cMap[5], linestyle='dashed')
            axs[1,2].set_xlabel('block ({n} trials)'.format(n=sSegments))
            axs[1,2].spines['top'].set_visible(False); axs[1,2].spines['right'].set_visible(False);
            axs[1,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
            axs[1,2].set_yticks(np.linspace(0, 1, 11)); 
                     
        fig.suptitle('From rev phase : {mouse}'.format(mouse=thisMouse))   
        plt.savefig('{directory}revPhaseAccuracyBinned_firstSix_{mouse}.pdf'.format(directory=outputPath, mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
        
    #accuracy segments from LAST 6 days of reversal phase
    if 4 in phaseIndex:
        if len(nRev)>6:
            nRev = np.round(np.linspace(phaseIndex.index(4),(len(phaseIndex) - 1 - phaseIndex[::-1].index(4)),(1+(len(phaseIndex) - 1 - phaseIndex[::-1].index(4))-phaseIndex.index(4)))).astype(int)
            cMap=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']       
            
            rateAccuracyGo_segment = [[] for _ in range(len(nRev))]; rateAccuracyNo_segment = [[] for _ in range(len(nRev))]
            for d in range(len(nRev)):
                    for s in range(len(iSegment[nRev[d]])):
                        if sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==3)>0:
                            rateAccuracyGo_segment[d].append(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)/(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==1)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==3)))
                        if sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==4)>0:
                            rateAccuracyNo_segment[d].append(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)/(sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==2)+sum(responseMatrix[nRev[d]][iSegment[nRev[d]][s]]==4)))
                                               
            fig,axs=plt.subplots(2,3,figsize=(10, 10))
            
            #
            if len(nRev)>=1:
                axs[0,0].plot(rateAccuracyGo_segment[len(nRev)-6], color=cMap[0], label=nRev[len(nRev)-6]+1)
                axs[0,0].plot(rateAccuracyNo_segment[len(nRev)-6], color=cMap[0], linestyle='dashed')
                axs[0,0].set_ylabel('accuracy (%)')
                axs[0,0].spines['top'].set_visible(False); axs[0,0].spines['right'].set_visible(False);
                axs[0,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[0,0].set_yticks(np.linspace(0, 1, 11)); 
        
            #
            if len(nRev)>=2:
                axs[0,1].plot(rateAccuracyGo_segment[len(nRev)-5], color=cMap[1], label=nRev[len(nRev)-5]+1)
                axs[0,1].plot(rateAccuracyNo_segment[len(nRev)-5], color=cMap[1], linestyle='dashed')
                axs[0,1].spines['top'].set_visible(False); axs[0,1].spines['right'].set_visible(False);
                axs[0,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[0,1].set_yticks(np.linspace(0, 1, 11)); 
        
            #
            if len(nRev)>=3:
                axs[0,2].plot(rateAccuracyGo_segment[len(nRev)-4], color=cMap[2], label=nRev[len(nRev)-4]+1)
                axs[0,2].plot(rateAccuracyNo_segment[len(nRev)-4], color=cMap[2], linestyle='dashed')
                axs[0,2].spines['top'].set_visible(False); axs[0,2].spines['right'].set_visible(False);
                axs[0,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[0,2].set_yticks(np.linspace(0, 1, 11)); 
        
            #
            if len(nRev)>=4:
                axs[1,0].plot(rateAccuracyGo_segment[len(nRev)-3], color=cMap[3], label=nRev[len(nRev)-3]+1)
                axs[1,0].plot(rateAccuracyNo_segment[len(nRev)-3], color=cMap[3], linestyle='dashed')
                axs[1,0].set_xlabel('block ({n} trials)'.format(n=sSegments))
                axs[1,0].set_ylabel('accuracy (%)')
                axs[1,0].spines['top'].set_visible(False); axs[1,0].spines['right'].set_visible(False);
                axs[1,0].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[1,0].set_yticks(np.linspace(0, 1, 11)); 
        
            #
            if len(nRev)>=5:
                axs[1,1].plot(rateAccuracyGo_segment[len(nRev)-2], color=cMap[4], label=nRev[len(nRev)-2]+1)
                axs[1,1].plot(rateAccuracyNo_segment[len(nRev)-2], color=cMap[4], linestyle='dashed')
                axs[1,1].set_xlabel('block ({n} trials)'.format(n=sSegments))
                axs[1,1].spines['top'].set_visible(False); axs[1,1].spines['right'].set_visible(False);
                axs[1,1].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[1,1].set_yticks(np.linspace(0, 1, 11)); 
        
            #
            if len(nRev)>=6:
                axs[1,2].plot(rateAccuracyGo_segment[len(nRev)-1], color=cMap[5], label=nRev[len(nRev)-1]+1)
                axs[1,2].plot(rateAccuracyNo_segment[len(nRev)-1], color=cMap[5], linestyle='dashed')
                axs[1,2].set_xlabel('block ({n} trials)'.format(n=sSegments))
                axs[1,2].spines['top'].set_visible(False); axs[1,2].spines['right'].set_visible(False);
                axs[1,2].legend(loc='upper left', bbox_to_anchor=(0, 1.1))
                axs[1,2].set_yticks(np.linspace(0, 1, 11)); 
                         
            fig.suptitle('From rev phase : {mouse}'.format(mouse=thisMouse))   
            plt.savefig('{directory}revPhaseAccuracyBinned_lastSix_{mouse}.pdf'.format(directory=outputPath, mouse=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            
                         
    
    ## binomial testing for response accuracy 
    # Ho == rateAccuracyGo (or rateAccuracyNo) = .50 (50%)
    # p<0.05
    """
    from scipy.stats import binomtest
    
    stats_go = []; stats_no = []; pVal_go = []; pVal_no = []; H_go = []; H_no = []; 
    for day in range(len(allDataSheet)):
        if phaseIndex[day] > 1:
            stats_go.append(binomtest(nHit[day], nHit[day]+nMiss[day], p=0.3, alternative='greater'))
            pVal_go.append(np.round(stats_go[day].pvalue,8))
            H_go.append(pVal_go[day]<0.05)
        elif phaseIndex[day] == 1:
            pVal_go.append('NaN')
            H_go.append('NaN')
            pVal_no.append('NaN')
            H_no.append('NaN')
        if phaseIndex[day] > 2:    
            stats_no.append(binomtest(nCorrectReject[day], nCorrectReject[day]+nFalseAlarm[day], p=0.5, alternative='greater'))
            pVal_no.append(np.round(stats_no[day].pvalue,8))
            H_no.append(pVal_no[day]<0.05)
    """
            
    # look at response latencies across sessions : this will first be an average of the day per type and then a blocked average
    respLatency_segment = [[[tResponse[d][iSegment[d][s]]] for s in range(len(iSegment[d]))] for d in range(len(allDataSheet))]
    
    
#%% Part 3 - Save out key data into a .py
    allDataSheet.to_csv('{directory}rawDataFile_{mouse}.csv'.format(directory=outputPath,mouse=thisMouse))
     
    def isVar(name):
        return name.startswith('t') and (len(name)>1 and name[1].isupper())
    saveVar = {name: value for name, value in globals().items() if isVar(name)}
    def isVar(name):
        return name.startswith('n') and (len(name)>1 and name[1].isupper())
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    def isVar(name):
        return name.startswith('i') and (len(name)>1 and name[1].isupper())
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    def isVar(name):
        return name.startswith('psth')
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    def isVar(name):
        return name.startswith('rate')
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    def isVar(name):
        return name.startswith('raw')
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    def isVar(name):
        return name.startswith('phase')
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    
    with open('{directory}suppDataFile.pkl'.format(directory=outputPath), 'wb') as file:
        pickle.dump(saveVar, file)


#%% Part 4 - Group-wide plotting and analysis
## going to plot the total number of days in each phase (this is going to be across mice)
elif not userInputs['oneMouse']:
    gOutputPath = 'S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\groupAnalysis\\'
    txtFiles = [file for file in os.listdir(gOutputPath) if file.startswith(tuple('theseMice'))]

    with open(gOutputPath + '\\' + txtFiles[0],'r') as dataLines:
        allMice = dataLines.readlines()
    for line in range(len(allMice)):
        allMice[line] = allMice[line][0:3]

# load in mouse, pull phase info and data for first and last of each phase out
    gPhase={'firstSpout': [], 'lastSpout': [], 'Spout': [], 
           'firstGo': [], 'lastGo': [], 'Go': [], 
           'firstNo': [], 'lastNo': [], 'No': [], 
           'firstReversal': [], 'lastReversal': [], 'Reversal': [],
           'psth_firstGo': [], 'psth_lastGo': [], 
           'psth_firstNo': [], 'psth_lastNo': [], 
           'psth_firstReversal': [], 'psth_lastReversal': [],
           }
    
    dataStruct = {}
    
    for thisMouse in range(len(allMice)):
        with open('{directory}suppDataFile.pkl'.format(directory=outputPath + '\\' + allMice[thisMouse] + '\\bx\\'), 'rb') as file:
            suppVar = pickle.load(file)
            dataStruct[thisMouse] = suppVar            
        
##        
    for thisMouse in range(len(allMice)):    
        if 1 in suppVar['phaseIndex']:
            gPhase['firstSpout'].append(dataStruct[thisMouse]['phaseIndex'].index(1)); gPhase['lastSpout'].append((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(1)))
            gPhase['Spout'].append(((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(1)))-dataStruct[thisMouse]['phaseIndex'].index(1))
        if 2 in suppVar['phaseIndex']:
            gPhase['firstGo'].append(dataStruct[thisMouse]['phaseIndex'].index(2)); gPhase['lastGo'].append((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(2)))
            gPhase['Go'].append(((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(2)))-dataStruct[thisMouse]['phaseIndex'].index(2))
        if 3 in suppVar['phaseIndex']:
            gPhase['firstNo'].append(dataStruct[thisMouse]['phaseIndex'].index(3)); gPhase['lastNo'].append((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3)))
            gPhase['No'].append(((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3)))-dataStruct[thisMouse]['phaseIndex'].index(3))
        if 4 in suppVar['phaseIndex']:
            gPhase['firstReversal'].append(dataStruct[thisMouse]['phaseIndex'].index(4)); gPhase['lastReversal'].append((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4)))
            gPhase['Reversal'].append(((len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4)))-dataStruct[thisMouse]['phaseIndex'].index(4))

## going to plot the accuracy from the first day of each phase across mice
    gFirstSesh=[]; gAllSeshAcc=[]; gFirstSeshAcc=[[] for mouse in range(len(allMice))]; groupAcc=[[] for day in range(1,4)];
    for thisMouse in range(len(allMice)):
        gFirstSesh.append(dataStruct[thisMouse]['iFirst'])
        gAllSeshAcc.append(dataStruct[thisMouse]['rateAccuracy_segment']) 
    for thisMouse in range(len(allMice)):
        for day in range(1,len(gFirstSesh[thisMouse])):
            gFirstSeshAcc[thisMouse].append(gAllSeshAcc[thisMouse][gFirstSesh[thisMouse][day]])
    for day in range(1,4):
        for thisMouse in range(len(allMice)):  
            if len(gFirstSesh[thisMouse])>day:
                groupAcc[day-1].append(gAllSeshAcc[thisMouse][gFirstSesh[thisMouse][day]][0])
            else:
                continue

    groupDayMean=[[] for day in range(1,4)]; groupDayStd=[[] for day in range(1,4)];        
    for day in range(len(groupAcc)):
        groupDayMean[day] = np.mean(groupAcc[day])
        groupDayStd[day] = np.std(groupAcc[day])    
        
    jitter=0.15; color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']       
    plt.figure(figsize=(14, 12))
    for day in range(len(groupAcc)):
        plt.errorbar(day+1,groupDayMean[day],yerr=groupDayStd[day],fmt=('_'),color='black',linewidth=3, capsize=10)
    for thisMouse in range(len(allMice)):
        for day in range(len(gFirstSeshAcc[thisMouse])):
            jitterAmt = day + 1 + np.random.uniform(-jitter,jitter)
            if day == 0:
                plt.scatter(jitterAmt,gFirstSeshAcc[thisMouse][day][0],c=color[thisMouse],s=100,label=allMice[thisMouse])
            else:
                plt.scatter(jitterAmt,gFirstSeshAcc[thisMouse][day][0],c=color[thisMouse],s=100)

    plt.legend(fontsize=24, title='mouse', loc='upper right', bbox_to_anchor=(1.1, 1));
    plt.xlabel('phase',fontsize=36)
    plt.ylabel('accuracy (%)',fontsize=36)
    plt.title('From first 20 trials of first session of each phase',fontsize=48)   
    ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    plt.tick_params(axis='both', labelsize=36); ax.set_xlim(0, 4); ax.set_ylim(0, 1); 
    ax.set_yticks(np.linspace(0, 1, 11)); ax.set_xticks([1,2,3],['go','no-go','reversal']); 
    plt.savefig('{directory}first20TrialAccuracy_firstSessionOfPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
           
    # group-wise phase accuracy

    gRateAccuracyGo=[[] for mice in range(len(allMice))]
    for thisMouse in range(len(allMice)):    
        nGo = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(2),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(2)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(2))-dataStruct[thisMouse]['phaseIndex'].index(2)))).astype(int)
        for i in range(len(nGo)):
            gRateAccuracyGo[thisMouse].append(dataStruct[thisMouse]['rateAccuracyGo'][nGo[i]])
    
    cMap=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']       
        
    plt.figure(figsize=(20, 16))
    plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
    plt.title('Accuracy across go sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
    plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks((range(1, max([[(len(gRateAccuracyGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
    ax.set_xlim(1, max([[(len(gRateAccuracyGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
    
    for thisMouse in range(len(allMice)):    
        plt.plot(list(range(1, len(gRateAccuracyGo[thisMouse])+1)), gRateAccuracyGo[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', label=allMice[thisMouse])
        
    plt.legend(fontsize=24, title='mouse', loc='upper right', bbox_to_anchor=(1.15, 1));
    
    plt.savefig('{directory}fullSessionAccuracy_goPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    ## no-go phase
    gRateAccuracyGo_noGo=[[] for mice in range(len(allMice))]
    gRateAccuracyNo_noGo=[[] for mice in range(len(allMice))]
    for thisMouse in range(len(allMice)):    
        if 3 in dataStruct[thisMouse]['phaseIndex']:
            nNo_noGo = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(3),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3))-dataStruct[thisMouse]['phaseIndex'].index(3)))).astype(int)
            nGo_noGo = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(3),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3))-dataStruct[thisMouse]['phaseIndex'].index(3)))).astype(int)
            for i in range(len(nNo_noGo)):
                gRateAccuracyGo_noGo[thisMouse].append(dataStruct[thisMouse]['rateAccuracyGo'][nGo_noGo[i]])
                gRateAccuracyNo_noGo[thisMouse].append(dataStruct[thisMouse]['rateAccuracyNo'][nNo_noGo[i]])
            
    # separate go and no go here (no-go dashed)      
    plt.figure(figsize=(20, 16))
    plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
    plt.title('Accuracy across no-go sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
    plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks((range(1, max([[(len(gRateAccuracyNo_noGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
    ax.set_xlim(1, max([[(len(gRateAccuracyNo_noGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
    
    for thisMouse in range(len(allMice)):    
        plt.plot(list(range(1, len(gRateAccuracyNo_noGo[thisMouse])+1)), gRateAccuracyNo_noGo[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', linestyle='dashed')
        plt.plot(list(range(1, len(gRateAccuracyGo_noGo[thisMouse])+1)), gRateAccuracyGo_noGo[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', label=allMice[thisMouse])
        
    plt.legend(fontsize=24,title='mouse', loc='upper right', bbox_to_anchor=(1.15, 1));
    
    plt.savefig('{directory}fullSessionAccuracy_noPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    
    ## reversal phase
    gRateAccuracyGo_rev=[[] for mice in range(len(allMice))]
    gRateAccuracyNo_rev=[[] for mice in range(len(allMice))]
    for thisMouse in range(len(allMice)):    
        if 4 in dataStruct[thisMouse]['phaseIndex']:
            nNo_rev = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(4),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4))-dataStruct[thisMouse]['phaseIndex'].index(4)))).astype(int)
            nGo_rev = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(4),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4))-dataStruct[thisMouse]['phaseIndex'].index(4)))).astype(int)
            for i in range(len(nNo_rev)):
                gRateAccuracyGo_rev[thisMouse].append(dataStruct[thisMouse]['rateAccuracyGo'][nGo_rev[i]])
                gRateAccuracyNo_rev[thisMouse].append(dataStruct[thisMouse]['rateAccuracyNo'][nNo_rev[i]])
            
    # separate go and no go here (no-go dashed)      
    plt.figure(figsize=(20, 16))
    plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
    plt.title('Accuracy across rev sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
    plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks((range(1, max([[(len(gRateAccuracyNo_rev[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
    ax.set_xlim(1, max([[(len(gRateAccuracyNo_rev[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
    
    for thisMouse in range(len(allMice)):    
        plt.plot(list(range(1, len(gRateAccuracyNo_rev[thisMouse])+1)), gRateAccuracyNo_rev[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', linestyle='dashed')
        plt.plot(list(range(1, len(gRateAccuracyGo_rev[thisMouse])+1)), gRateAccuracyGo_rev[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', label=allMice[thisMouse])
        
    plt.legend(fontsize=24, title='mouse', loc='upper right', bbox_to_anchor=(1.15, 1));
    
    plt.savefig('{directory}fullSessionAccuracy_revPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    ## change in too early during go/no-go
    gRateTooEarly_noGo=[[] for mice in range(len(allMice))]
    for thisMouse in range(len(allMice)):    
        if 3 in dataStruct[thisMouse]['phaseIndex']:
            nNo_noGo = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(3),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(3))-dataStruct[thisMouse]['phaseIndex'].index(3)))).astype(int)
            for i in range(len(nNo_noGo)):
                gRateTooEarly_noGo[thisMouse].append(dataStruct[thisMouse]['rateTooEarly'][nNo_noGo[i]])
            
    # separate go and no go here (no-go dashed)      
    plt.figure(figsize=(20, 16))
    plt.title('TooEarly across no-go sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('TooEarly (%)',fontsize=36)
    plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks((range(1, max([[(len(gRateTooEarly_noGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,.6,7)); 
    ax.set_xlim(1, max([[(len(gRateTooEarly_noGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, .6)
    
    for thisMouse in range(len(allMice)):    
        plt.plot(list(range(1, len(gRateTooEarly_noGo[thisMouse])+1)), gRateTooEarly_noGo[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', label=allMice[thisMouse])
        
    plt.legend(fontsize=24,title='mouse', loc='upper right', bbox_to_anchor=(1.15, 1));
    
    plt.savefig('{directory}fullSessionTooEarly_noPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    ## change in too early during reversal
    gRateTooEarly_rev=[[] for mice in range(len(allMice))]
    for thisMouse in range(len(allMice)):    
        if 4 in dataStruct[thisMouse]['phaseIndex']:
            nRev = np.round(np.linspace(dataStruct[thisMouse]['phaseIndex'].index(4),(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4)),(1+(len(dataStruct[thisMouse]['phaseIndex']) - 1 - dataStruct[thisMouse]['phaseIndex'][::-1].index(4))-dataStruct[thisMouse]['phaseIndex'].index(4)))).astype(int)
            for i in range(len(nRev)):
                gRateTooEarly_rev[thisMouse].append(dataStruct[thisMouse]['rateTooEarly'][nRev[i]])
              
    plt.figure(figsize=(20, 16))
    plt.title('TooEarly across rev sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('TooEarly (%)',fontsize=36)
    plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
    ax.set_xticks((range(1, max([[(len(gRateTooEarly_rev[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,.6,7)); 
    ax.set_xlim(1, max([[(len(gRateTooEarly_rev[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, .6)
    
    for thisMouse in range(len(allMice)):    
        plt.plot(list(range(1, len(gRateTooEarly_rev[thisMouse])+1)), gRateTooEarly_rev[thisMouse], color=cMap[thisMouse], linewidth=4, markersize=10, marker='o', label=allMice[thisMouse])
        
    plt.legend(fontsize=24,title='mouse', loc='upper right', bbox_to_anchor=(1.15, 1));
    
    plt.savefig('{directory}fullSessionTooEarly_revPhase_allMice.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
    plt.show(); plt.close()
    
    ##
    def isVar(name):
        return name.startswith('g') and (len(name)>1 and name[1].isupper())
    saveVar = {name: value for name, value in globals().items() if isVar(name)}
    def isVar(name):
        return name.startswith('n') and (len(name)>1 and name[1].isupper())
    saveVar.update({name: value for name, value in globals().items() if isVar(name)})
    
    with open('{directory}groupDataFile.pkl'.format(directory=gOutputPath), 'wb') as file:
        pickle.dump(saveVar, file)
        
#%% Part 5 - analysis post-unblinding (for groups/averages)       
    ### do group analysis after unblinding (put txt doc titled 'doGroups' in groupAnalysis folder)
    
    #check doGroups.txt is in the folder of interest
    gOutputPath = 'S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\groupAnalysis\\'
    txtFiles = [file for file in os.listdir(gOutputPath) if file.startswith(tuple('doGroups'))]
    if txtFiles:
        with open(gOutputPath + '\\' + txtFiles[0],'r') as dataLines:
            theseGroups = dataLines.readlines()
        groupedMice = []
        for line in theseGroups:
            condition, mice_str = line.strip().split(':')
            mice = mice_str.split(',')
            groupedMice.append([condition]+mice)
    
        #group initial accuracy from each phase
        groupDayMean=[[[] for day in range(1,4)] for group in range(len(groupedMice))]; groupDayStd=[[[] for day in range(1,4)] for group in range(len(groupedMice))];        
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            groupIDs = []
            for thisMouse in range(len(theseMice)):
                groupIDs.append(allMice.index(theseMice[thisMouse]))
                for day in range(len(groupAcc)):
                    groupDayMean[thisGroup][day] = np.mean([groupAcc[day][i] for i in groupIDs])
                    groupDayStd[thisGroup][day] = np.std([groupAcc[day][i] for i in groupIDs])    
                    
        cMap=['#a9a9a9', '#8b008b', '#17becf', '#bcbd22']; conditions=['go','no-go','reversal']   
        plt.figure(figsize=(14, 12))
        barWidth = 0.45
        widthBtwDays = .01
        r1 = [x + widthBtwDays for x in np.arange(len(groupDayMean[0]))]
        r2 = [x + barWidth for x in r1]
        
        # Plotting bars
        plt.bar(r1, groupDayMean[0], color=cMap[0], width=barWidth, capsize=7.5, label=groupedMice[0][0])
        plt.bar(r2, groupDayMean[1], color=cMap[1], width=barWidth, capsize=7.5, label=groupedMice[1][0])
        
        # Adding error bars separately with customized elinewidth
        plt.errorbar(r1, groupDayMean[0], yerr=groupDayStd[0], fmt='none', capsize=7.5, elinewidth=3, color='black')
        plt.errorbar(r2, groupDayMean[1], yerr=groupDayStd[1], fmt='none', capsize=7.5, elinewidth=3, color='black')

        plt.xticks([r + widthBtwDays + barWidth/2 for r in range(len(groupDayMean[0]))], conditions)
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.1, 1));
        plt.xlabel('phase',fontsize=36)
        plt.ylabel('accuracy (%)',fontsize=36)
        plt.title('From start of each phase',fontsize=36)   
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_yticks(np.linspace(0, 1, 11));  
        plt.tick_params(axis='both', labelsize=36); ax.set_ylim(0, 0.7); 
        plt.savefig('{directory}first20TrialAccuracy_firstSessionOfPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
         
        #grouped go phase accuracy plot
        groupedRateAccuracyGo=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            for thisMouse in range(len(theseMice)):
                groupedRateAccuracyGo[thisGroup][thisMouse] = gRateAccuracyGo[allMice.index(theseMice[thisMouse])]
                
        plt.figure(figsize=(20, 16))
        plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
        plt.title('Accuracy across go sessions',fontsize=36); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
        plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xticks((range(1, max([[(len(gRateAccuracyGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
        ax.set_xlim(1, max([[(len(gRateAccuracyGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
        for thisGroup in range(len(groupedMice)):
            for theseMice in range(len(groupedMice[thisGroup][1:])):
                mouseID = groupedMice[thisGroup][1+theseMice]
                if theseMice==0:
                    plt.plot(list(range(1, len(gRateAccuracyGo[allMice.index(mouseID)])+1)), gRateAccuracyGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                else:
                    plt.plot(list(range(1, len(gRateAccuracyGo[allMice.index(mouseID)])+1)), gRateAccuracyGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o',)
                
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.15, 1));
        plt.savefig('{directory}fullSessionAccuracy_goPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.show(); plt.close()
        
        ##grouped no-go phase accuracy plot
        groupedRateAccuracyGo_noGo=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        groupedRateAccuracyNo_noGo=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            for thisMouse in range(len(theseMice)):
                groupedRateAccuracyGo_noGo[thisGroup][thisMouse] = gRateAccuracyGo_noGo[allMice.index(theseMice[thisMouse])]
                groupedRateAccuracyNo_noGo[thisGroup][thisMouse] = gRateAccuracyNo_noGo[allMice.index(theseMice[thisMouse])]
                
        plt.figure(figsize=(20, 16))
        plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
        plt.title('Accuracy across no-go sessions',fontsize=36); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
        plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xticks((range(1, max([[(len(gRateAccuracyNo_noGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
        ax.set_xlim(1, max([[(len(gRateAccuracyNo_noGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
        for thisGroup in range(len(groupedMice)):
            for theseMice in range(len(groupedMice[thisGroup][1:])):
                mouseID = groupedMice[thisGroup][1+theseMice]
                if theseMice==0:
                    plt.plot(list(range(1, len(gRateAccuracyGo_noGo[allMice.index(mouseID)])+1)), gRateAccuracyGo_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                    plt.plot(list(range(1, len(gRateAccuracyNo_noGo[allMice.index(mouseID)])+1)), gRateAccuracyNo_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linestyle='dashed', linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                else:
                    plt.plot(list(range(1, len(gRateAccuracyGo_noGo[allMice.index(mouseID)])+1)), gRateAccuracyGo_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o',)
                    plt.plot(list(range(1, len(gRateAccuracyNo_noGo[allMice.index(mouseID)])+1)), gRateAccuracyNo_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linestyle='dashed', linewidth=4, markersize=10, marker='o',)
                
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.20, 1));
        plt.savefig('{directory}fullSessionAccuracy_noGoPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.show(); plt.close()
        
        ##grouped reversal phase accuracy plot
        groupedRateAccuracyGo_rev=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        groupedRateAccuracyNo_rev=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            for thisMouse in range(len(theseMice)):
                groupedRateAccuracyGo_rev[thisGroup][thisMouse] = gRateAccuracyGo_rev[allMice.index(theseMice[thisMouse])]
                groupedRateAccuracyNo_rev[thisGroup][thisMouse] = gRateAccuracyNo_rev[allMice.index(theseMice[thisMouse])]
                
        plt.figure(figsize=(20, 16))
        plt.axhline(0.65, color='black', linestyle='--', linewidth=2)
        plt.title('Accuracy across no-go sessions',fontsize=36); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('Accuracy (%)',fontsize=36)
        plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xticks((range(1, max([[(len(gRateAccuracyNo_rev[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,1,11)); 
        ax.set_xlim(1, max([[(len(gRateAccuracyNo_rev[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, 1)
        for thisGroup in range(len(groupedMice)):
            for theseMice in range(len(groupedMice[thisGroup][1:])):
                mouseID = groupedMice[thisGroup][1+theseMice]
                if theseMice==0:
                    plt.plot(list(range(1, len(gRateAccuracyGo_rev[allMice.index(mouseID)])+1)), gRateAccuracyGo_rev[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                    plt.plot(list(range(1, len(gRateAccuracyNo_rev[allMice.index(mouseID)])+1)), gRateAccuracyNo_rev[allMice.index(mouseID)], color=cMap[thisGroup], linestyle='dashed', linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                else:
                    plt.plot(list(range(1, len(gRateAccuracyGo_rev[allMice.index(mouseID)])+1)), gRateAccuracyGo_rev[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o',)
                    plt.plot(list(range(1, len(gRateAccuracyNo_rev[allMice.index(mouseID)])+1)), gRateAccuracyNo_rev[allMice.index(mouseID)], color=cMap[thisGroup], linestyle='dashed', linewidth=4, markersize=10, marker='o',)
                
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.20, 1));
        plt.savefig('{directory}fullSessionAccuracy_revPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.show(); plt.close()
        
        ## change in too early during go/no-go
        groupedRateTooEarly_noGo=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            for thisMouse in range(len(theseMice)):
                groupedRateTooEarly_noGo[thisGroup][thisMouse] = gRateTooEarly_noGo[allMice.index(theseMice[thisMouse])]
                
        # separate go and no go here (no-go dashed)      
        plt.figure(figsize=(20, 16))
        plt.title('TooEarly across no-go sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('TooEarly (%)',fontsize=36)
        plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xticks((range(1, max([[(len(gRateTooEarly_noGo[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,.6,7)); 
        ax.set_xlim(1, max([[(len(gRateTooEarly_noGo[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, .6)
        
        for thisGroup in range(len(groupedMice)):
            for theseMice in range(len(groupedMice[thisGroup][1:])):
                mouseID = groupedMice[thisGroup][1+theseMice]
                if theseMice==0:
                    plt.plot(list(range(1, len(gRateTooEarly_noGo[allMice.index(mouseID)])+1)), gRateTooEarly_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                else:
                    plt.plot(list(range(1, len(gRateTooEarly_noGo[allMice.index(mouseID)])+1)), gRateTooEarly_noGo[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o',)
         
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.2, 1));
        
        plt.savefig('{directory}fullSessionTooEarly_noPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.show(); plt.close()
        
        ## change in too early during reversal
        groupedRateTooEarly_rev=[[[] for group in range(len(groupedMice))] for mice in range(len(groupedMice[0][1:]))]
        for thisGroup in range(len(groupedMice)):
            theseMice = groupedMice[thisGroup][1:]
            for thisMouse in range(len(theseMice)):
                groupedRateTooEarly_rev[thisGroup][thisMouse] = gRateTooEarly_rev[allMice.index(theseMice[thisMouse])]
                
        # separate go and no go here (no-go dashed)      
        plt.figure(figsize=(20, 16))
        plt.title('TooEarly across no-go sessions',fontsize=48); plt.xlabel('Session (#)',fontsize=36); plt.ylabel('TooEarly (%)',fontsize=36)
        plt.tick_params(axis='both', labelsize=36); ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_xticks((range(1, max([[(len(gRateTooEarly_rev[m]))] for m in range(len(allMice))])[0]+1))); ax.set_yticks(np.linspace(0,.6,7)); 
        ax.set_xlim(1, max([[(len(gRateTooEarly_rev[m]))] for m in range(len(allMice))])[0]); ax.set_ylim(0, .6)
        
        for thisGroup in range(len(groupedMice)):
            for theseMice in range(len(groupedMice[thisGroup][1:])):
                mouseID = groupedMice[thisGroup][1+theseMice]
                if theseMice==0:
                    plt.plot(list(range(1, len(gRateTooEarly_rev[allMice.index(mouseID)])+1)), gRateTooEarly_rev[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o', label=groupedMice[thisGroup][0])
                else:
                    plt.plot(list(range(1, len(gRateTooEarly_rev[allMice.index(mouseID)])+1)), gRateTooEarly_rev[allMice.index(mouseID)], color=cMap[thisGroup], linewidth=4, markersize=10, marker='o',)
         
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.2, 1));
        
        plt.savefig('{directory}fullSessionTooEarly_noPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.show(); plt.close()
 

        #group first-last day accuracy in each phase
        groupedFirstLastSessionAcc = [[] for group in range(len(groupedMice))]; groupedFirstLastSessionAccMean = [[] for group in range(len(groupedMice))]; groupedFirstLastSessionAccStd = [[] for group in range(len(groupedMice))];
        for thisGroup in range(len(groupedMice)):
            groupedFirstLastSessionAcc[thisGroup]=[[[] for mouse in range(len(groupedMice[thisGroup][1:]))] for day in range(1,7)];
            groupedFirstLastSessionAccMean[thisGroup]=[[] for day in range(1,7)]; 
            groupedFirstLastSessionAccStd[thisGroup]=[[] for day in range(1,7)];
            theseMice = groupedMice[thisGroup][1:]
            groupIDs = []; groupDates = [[] for mouse in range(len(groupedMice))];
            for thisMouse in range(len(theseMice)):
                groupIDs.append(allMice.index(theseMice[thisMouse]))
                groupDates[thisMouse] = [gPhase['firstGo'][groupIDs[thisMouse]],gPhase['lastGo'][groupIDs[thisMouse]],gPhase['firstNo'][groupIDs[thisMouse]],gPhase['lastNo'][groupIDs[thisMouse]],gPhase['firstReversal'][groupIDs[thisMouse]],gPhase['lastReversal'][groupIDs[thisMouse]]]
                for day in range(len(groupedFirstLastSessionAccMean[thisGroup])):
                    groupedFirstLastSessionAcc[thisGroup][day][thisMouse] = np.mean(gAllSeshAcc[groupIDs[thisMouse]][groupDates[thisMouse][day]])
    
            for day in range(len(groupedFirstLastSessionAccMean[thisGroup])):
                groupedFirstLastSessionAccMean[thisGroup][day] = np.mean(groupedFirstLastSessionAcc[thisGroup][day])
                groupedFirstLastSessionAccStd[thisGroup][day] = np.std(groupedFirstLastSessionAcc[thisGroup][day])
                
        cMap=['#a9a9a9', '#8b008b', '#17becf', '#bcbd22']; conditions=['first go','last go','first no-go','last no-go','first reversal','last reversal']   
        plt.figure(figsize=(14, 12))
        barWidth = 0.45
        widthBtwSessions = .01
        r1 = [x + widthBtwSessions for x in np.arange(len(groupedFirstLastSessionAccMean[0]))]
        r2 = [x + barWidth for x in r1]
        
        # Plotting bars
        plt.bar(r1, groupedFirstLastSessionAccMean[0], color=cMap[0], width=barWidth, capsize=7.5, label=groupedMice[0][0])
        plt.bar(r2, groupedFirstLastSessionAccMean[1], color=cMap[1], width=barWidth, capsize=7.5, label=groupedMice[1][0])
        
        # Adding error bars separately with customized elinewidth
        plt.errorbar(r1, groupedFirstLastSessionAccMean[0], yerr=groupedFirstLastSessionAccStd[0], fmt='none', capsize=7.5, elinewidth=3, color='black')
        plt.errorbar(r2, groupedFirstLastSessionAccMean[1], yerr=groupedFirstLastSessionAccStd[1], fmt='none', capsize=7.5, elinewidth=3, color='black')

        plt.xticks([r + widthBtwSessions + barWidth/2 for r in range(len(groupedFirstLastSessionAccMean[0]))], conditions)
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.1, 1));
        plt.xlabel('phase',fontsize=36)
        plt.ylabel('accuracy (%)',fontsize=36)
        plt.title('From start of each phase',fontsize=36)   
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        ax.set_yticks(np.linspace(0, 1, 11));  
        plt.tick_params(axis='both', labelsize=36); ax.set_ylim(0, 1); 
        plt.savefig('{directory}trialAccuracy_firstLastSessionsOfPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        
        #group first-last day accuracy in each phase
        groupedTimeToCriterion = [[] for group in range(len(groupedMice))]; indexFinalSession = [1,3,5];
        for thisGroup in range(len(groupedMice)):
            groupedTimeToCriterion[thisGroup]=[[] for mouse in range(len(groupedMice[thisGroup][1:]))];
            
            theseMice = groupedMice[thisGroup][1:]
            groupIDs = []; groupDates = [[] for mouse in range(len(groupedMice))];
            for thisMouse in range(len(theseMice)):
                groupIDs.append(allMice.index(theseMice[thisMouse]))
                groupDates[thisMouse] = [gPhase['firstGo'][groupIDs[thisMouse]],gPhase['lastGo'][groupIDs[thisMouse]],gPhase['firstNo'][groupIDs[thisMouse]],gPhase['lastNo'][groupIDs[thisMouse]],gPhase['firstReversal'][groupIDs[thisMouse]],gPhase['lastReversal'][groupIDs[thisMouse]]]
                
                for thisPhase in range(len(indexFinalSession)):
                    groupedTimeToCriterion[thisGroup][thisMouse].append(groupDates[thisMouse][indexFinalSession[thisPhase]]-groupDates[thisMouse][indexFinalSession[thisPhase]-1])
        
        groupedTimeToCriterion_mean = [[] for group in range(len(groupedMice))] 
        groupedTimeToCriterion_std = [[] for group in range(len(groupedMice))]
        for thisGroup in range(len(groupedMice)):
            for thisPhase in range(len(indexFinalSession)):
                tempArray = []
                for thisMouse in range(len(groupedMice[thisGroup][1::])):
                    tempArray.append(groupedTimeToCriterion[thisGroup][thisMouse][thisPhase])
                    
                groupedTimeToCriterion_mean[thisGroup].append(np.mean(tempArray))
                groupedTimeToCriterion_std[thisGroup].append(np.std(tempArray)) 
        

        cMap=['#a9a9a9', '#8b008b', '#17becf', '#bcbd22']; conditions=['go','no-go','reversal']   
        plt.figure(figsize=(14, 12))
        barWidth = 0.45
        widthBtwSessions = .01
        r1 = [x + widthBtwSessions for x in np.arange(len(groupedTimeToCriterion_mean[0]))]
        r2 = [x + barWidth for x in r1]
        
        # Plotting bars
        plt.bar(r1, groupedTimeToCriterion_mean[0], color=cMap[0], width=barWidth, capsize=7.5, label=groupedMice[0][0])
        plt.bar(r2, groupedTimeToCriterion_mean[1], color=cMap[1], width=barWidth, capsize=7.5, label=groupedMice[1][0])
        
        # Adding error bars separately with customized elinewidth
        plt.errorbar(r1, groupedTimeToCriterion_mean[0], yerr=groupedTimeToCriterion_std[0], fmt='none', capsize=7.5, elinewidth=3, color='black')
        plt.errorbar(r2, groupedTimeToCriterion_mean[1], yerr=groupedTimeToCriterion_std[1], fmt='none', capsize=7.5, elinewidth=3, color='black')

        plt.xticks([r + widthBtwSessions + barWidth/2 for r in range(len(groupedTimeToCriterion_mean[0]))], conditions)
        plt.legend(fontsize=24, loc='upper right', bbox_to_anchor=(1.1, 1));
        plt.xlabel('phase',fontsize=36)
        plt.ylabel('trials to criteria (#)',fontsize=36)
        plt.title('Number of sessions to criteria',fontsize=36)   
        ax = plt.gca(); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False);
        plt.tick_params(axis='both', labelsize=36); 
        plt.savefig('{directory}trialAccuracy_firstLastSessionsOfPhase_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
        plt.savefig('{directory}numberOfTrialsToCriteria_grouped.pdf'.format(directory=gOutputPath), format='pdf', bbox_inches='tight', dpi=300)
