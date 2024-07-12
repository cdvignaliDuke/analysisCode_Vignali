# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 12:35:25 2023

@author: Carlo Vignali

@def: *localFunctions* - group of basic local functions useful for my scripts
    
"""

#import packages and system modules
import shutil
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

#%% identifyOutliers
def identifyOutliers(arr, threshold):
    first_derivative = np.diff(arr, prepend=arr[0])  # Calculate the first derivative
    mean = np.mean(first_derivative)
    std_dev = np.std(first_derivative)
    outliers = np.abs(first_derivative - mean) > threshold * std_dev # Identify outliers
    return outliers

#%% findIndices
def findOutlierIndices(data, threshold):
    outlier_indices = []
    for i, trial in enumerate(data):
        outliers = identifyOutliers(trial, threshold)
        if np.any(outliers):
            outlier_indices.append(i)
    return outlier_indices

#%% replaceOutliers
def replaceOutliers(arr, outliers):
    arr = np.array(arr)  # Ensure input is a NumPy array

    for i in range(len(arr)):
        if outliers[i]:
            if i == 0:  # First element
                arr[i] = arr[i+1]
            elif i == len(arr) - 1:  # Last element
                arr[i] = arr[i-1]
            else:
                arr[i] = (arr[i-1] + arr[i+1]) / 2  # Replace outlier with average of neighbors

    return arr
    
#%% moveFiles - moves files from dump folder (mouse) into holding folder (yymmdd) for 
def moveFiles(thisMouse,inputPath,prefix='',suffix='',fileType='',outputPath=None): 
    date=[]
    if outputPath is None: outputPath = inputPath + thisMouse + '\\'
    toMove = [file for file in os.listdir(inputPath) if '{p}{m}{s}'.format(m=thisMouse,p=prefix,s=suffix) in file and file.endswith(fileType)] #and os.path.isfile(os.path.join(inputPath, file))]
    if not toMove == []:
        for file in range(len(toMove)):
            if toMove[file].endswith('.txt'):
                date = toMove[file][2:4]+toMove[file][5:7]+toMove[file][8:10]
                shutil.move(inputPath + toMove[file],outputPath + date)
            elif toMove[file].endswith('.avi'):
                shutil.move(inputPath + toMove[file],outputPath)
            elif toMove[file].endswith('.mp4'):
                date = toMove[file][12:14]+toMove[file][15:17]+toMove[file][18:20]
                shutil.move(inputPath + toMove[file],outputPath + date)
            elif os.path.isdir(os.path.join(inputPath,toMove[file])):
                date = toMove[file][4:10]
                shutil.move(inputPath + toMove[file],outputPath + date)  
                
    
#%% getUserInput - this function is designed to receive an input question (prompt), look for a boolean response from the user (key), 
# and use that response to determine decision forks in larger scripts

def getUserInput(prompt, key):        
    response = {}
    while True:
        try:
            value = input(prompt)
            response[key] = bool(int(value))
            break
        except ValueError:
            print('Invalid input. Please enter 1|0.')
 
    try:
        userInputs[key] = response[key]
    except NameError:
        userInputs = {}
        userInputs[key] = response[key]
       
    return userInputs    

#%% exportMPC - this (beautiful) function is designed to import the MPC text file, separate the variables into columns to be fit into an excel 
#and then provide a variable list associated with the columns for that excel [as of right now that list looks like @line 19]
    
def exportMPC(dataPath):  

    dataSheet = pd.DataFrame(columns = ['Subject', 'Experiment', 'Date', 'startTime', 'endTime', 'totalTime', 'ITIRange', 'maxTrials', 'endTrial', 'trialStartTime', 'trialEndTime', 'trialType', 'trialITI', 'responseTime', 'rewardTime', 'responseCounter', 'toneTime', 'terminalBeamBreakTime', 'beamBreakTime', 'trialOutcome'])
    
    txtFiles = [file for file in os.listdir(dataPath) if file.endswith(tuple('.txt'))]
    movFiles = [file for file in os.listdir(dataPath) if file.endswith(tuple('.mp4'))]
    if not txtFiles and movFiles:
        print('Script closing. No txt file found in folder: {}'.format(dataPath))
        sys.exit()
    elif not txtFiles and not movFiles:
        return 

    with open(dataPath + '\\' + txtFiles[0],'r') as dataLines:
        MPCArray = dataLines.readlines()
        MPCArray = MPCArray[1:]
        
        # !! ADD DATE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'Start Date:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'Date'] = (MPCArray[itemIndex][-3:-1]+MPCArray[itemIndex][-9:-7]+MPCArray[itemIndex][-6:-4])
        else:
            print('No "Date" found for file. Program will continue, but make note that file may be improperly formatted.')    
            
        # !! ADD SUBJECT NUMBER AND EXPERIMENT !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'Subject:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'Subject'] = (MPCArray[itemIndex][-4:-1])
            dataSheet.at[0,'Experiment'] = (MPCArray[itemIndex][-12:-5])
                        
        # !! ADD START TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'Start Time:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'startTime'] = (MPCArray[itemIndex][-9:-1])
            
        # !! ADD END TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'End Time:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'endTime'] = (MPCArray[itemIndex][-9:-1]) 
            
        # !! ADD TOTAL TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'T:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for dot in range(len(dotIndex)):
                if tempData[dotIndex[dot]-2]==' ':
                    dotArray.append(tempData[dotIndex[dot]-1:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-3]==' ':
                    dotArray.append(tempData[dotIndex[dot]-2:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-4]==' ':
                    dotArray.append(tempData[dotIndex[dot]-3:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-5]==' ':
                    dotArray.append(tempData[dotIndex[dot]-4:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'totalTime'] = np.array(altDotArray)     
            
        # !! ADD MAX TRIALS !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'N:  ' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'maxTrials'] = np.array(float(MPCArray[itemIndex][-8:-1]))

        # !! ADD TOTAL TRIALS !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'W:  ' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'endTrial'] = np.array(float(MPCArray[itemIndex][-8:-1]))  
            
        # !! ADD ITI RANGE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'Y:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1]
            dot = '.'
            dotIndex = []
            for index, char in enumerate(tempData):
                if char == dot:
                    dotIndex.append(index)

            dotArray = []
            for dot in range(len(dotIndex)):
                if tempData[dotIndex[dot]-3]==' ':
                    dotArray.append(tempData[dotIndex[dot]-2:dotIndex[dot]+4])
                elif not tempData[dotIndex[dot]-3]==' ':
                    dotArray.append(tempData[dotIndex[dot]-3:dotIndex[dot]+4])

            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'ITIRange'] = np.array(altDotArray)
            
        # !! ADD TRIAL START TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'A: ' in row:
                break
            elif 'A:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialStartTime'] = np.array(altDotArray)         

        # !! ADD TRIAL END TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'B: ' in row:
                break
            elif 'B:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialEndTime'] = np.array(altDotArray)         

        # !! ADD TRIAL TYPE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'C: ' in row:
                break
            elif 'C:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialType'] = np.array(altDotArray)         

        # !! ADD TRIAL ITI !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'D: ' in row:
                break
            elif 'D:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialITI'] = np.array(altDotArray)         

        # !! ADD TRIAL RESPONSE TIME CLOSE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'E: ' in row:
                break
            elif 'E:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'responseTime'] = np.array(altDotArray)         

        # !! ADD TERMINAL BEAM BREAK - THESE ARE SPECIFICALLY BEAM BREAKS THAT ENDED THE TRIAL !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'F: ' in row:
                break
            elif 'F:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'terminalBeamBreakTime'] = np.array(altDotArray)         

        # !! ADD REWARD TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'G: ' in row:
                break
            elif 'G:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'rewardTime'] = np.array(altDotArray)         

        # !! ADD RESPONSE COUNTER !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'I: ' in row:
                break
            elif 'I:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1]
            dot = '.'
            dotIndex = []
            for index, char in enumerate(tempData):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for dot in range(len(dotIndex)):
                if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData):
                    if tempData[dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[dotIndex[dot]-1:dotIndex[dot]+4])
                    elif tempData[dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[dotIndex[dot]-4:dotIndex[dot]+4])
                else:
                    break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'responseCounter'] = np.array(altDotArray)         

        # !! ADD TRIAL OUTCOME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'K: ' in row:
                break
            elif 'K:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialOutcome'] = np.array(altDotArray)         

        # !! ADD TRIAL OUTCOME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'J: ' in row:
                break
            elif 'J:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+(int(math.ceil(int(dataSheet.loc[0,'endTrial'])/5)))+1]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'toneTime'] = np.array(altDotArray)         

        # !! ADD BEAM BREAKS !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'M: ' in row:
                break
            elif 'M:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+1000]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if dotIndex[dot] - 3 >= 0 and dotIndex[dot] + 4 <= len(tempData[row]):
                        if tempData[row][dotIndex[dot]-2]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-3]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-4]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                        elif tempData[row][dotIndex[dot]-5]==' ':
                            dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    else:
                        break
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'beamBreakTime'] = np.array(altDotArray)        

    return dataSheet                    
    dataSheet.to_excel('{}\\{}.xlsx'.format(dataPath,txtFiles[0][:-4]), index=True)

