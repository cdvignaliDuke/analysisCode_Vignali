# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 13:56:26 2023

@author: Carlo Vignali

@def: *exportMPC*
    - this function is designed to import the MPC text file, separate the variables into columns to be fit into an excel 
and then provide a variable list associated with the columns for that excel [as of right now that list looks like @line 19]
    
"""
    'Subject' : [],
    'Experiment' : [],
    'Group' : [],
    'Date' : [],
    'maxTrials' : [],
    'ITIRange' : [],
    'endTime' : [],
    'endTrial' : [],
    'trialStartTime' : [],
    'trialEndTime' : [],
    'trialType' : [],
    'trialITI' : [],
    'responseTime' : [],
    'rewardTime' : [],
    'responseCounter' : [],
    'toneTime' : [],
    'trialResult' : [],
    'beamBreakTime' : [] 
    
    }

import os
import sys


#    
def exportMPC(dataPath):  

    dataSheet = pd.DataFrame(columns = ['Subject', 'Experiment', 'Date', 'startTime', 'endTime', 'totalTime', 'ITIRange', 'maxTrials', 'endTrial', 'trialStartTime', 'trialEndTime', 'trialType', 'trialITI', 'responseTime', 'rewardTime', 'responseCounter', 'toneTime', 'terminalBeamBreakTime', 'beamBreakTime', 'trialOutcome'])
    
    currFolder = os.listdir(dataPath)
    txtFiles = [file for file in os.listdir(dataPath) if file.endswith(tuple('.txt'))]
    if not txtFiles:
        print('Script closing. No txt file found in folder: {}'.format(dataPath))
        sys.exit()

    with open(dataPath + '\\' + currFolder[0],'r') as dataLines:
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
            dataSheet.at[0,'Date'] = np.array(MPCArray[itemIndex][-9:])   
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
            dataSheet.at[0,'Subject'] = np.array(MPCArray[itemIndex][-4:-1])
            dataSheet.at[0,'Experiment'] = np.array(MPCArray[itemIndex][-12:-5])
                        
        # !! ADD START TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'Start Time:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'startTime'] = np.array(MPCArray[itemIndex][-9:-1])
            
        # !! ADD END TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'End Time:' in row:
                itemIndex = line
                found = True
                break
        if found:
            dataSheet.at[0,'endTime'] = np.array(MPCArray[itemIndex][-9:-1]) 
            
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
            if tempData[dotIndex[0]-5]==' ':
                dotArray.append(tempData[dotIndex[0]-4:dotIndex[0]+4])
            elif tempData[dotIndex[dot]-4]==' ':
                dotArray.append(tempData[dotIndex[0]-3:dotIndex[0]+4])
            elif tempData[dotIndex[dot]-3]==' ':
                dotArray.append(tempData[dotIndex[0]-2:dotIndex[0]+4])
            elif tempData[dotIndex[dot]-2]==' ':
                dotArray.append(tempData[dotIndex[0]-1:dotIndex[0]+4])
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
            if 'A:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialStartTime'] = np.array(altDotArray)         

        # !! ADD TRIAL END TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'B:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialEndTime'] = np.array(altDotArray)         

        # !! ADD TRIAL TYPE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'C:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialType'] = np.array(altDotArray)         

        # !! ADD TRIAL ITI !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'D:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialITI'] = np.array(altDotArray)         

        # !! ADD TRIAL RESPONSE TIME CLOSE !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'E:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'responseTime'] = np.array(altDotArray)         

        # !! ADD TERMINAL BEAM BREAK - THESE ARE SPECIFICALLY BEAM BREAKS THAT ENDED THE TRIAL !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'F:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'terminalBeamBreakTime'] = np.array(altDotArray)         

        # !! ADD REWARD TIME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'G:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'rewardTime'] = np.array(altDotArray)         

        # !! ADD RESPONSE COUNTER !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'I:' in row:
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
                if tempData[dotIndex[dot]-5]==' ':
                    dotArray.append(tempData[dotIndex[dot]-4:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-4]==' ':
                    dotArray.append(tempData[dotIndex[dot]-3:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-3]==' ':
                    dotArray.append(tempData[dotIndex[dot]-2:dotIndex[dot]+4])
                elif tempData[dotIndex[dot]-2]==' ':
                    dotArray.append(tempData[dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'responseCounter'] = np.array(altDotArray)         

        # !! ADD TRIAL OUTCOME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'K:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'trialOutcome'] = np.array(altDotArray)         

        # !! ADD TRIAL OUTCOME !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'J:' in row:
                itemIndex = line
                found = True
                break
        if found:
            tempData = MPCArray[itemIndex+1:itemIndex+37]
            dot = '.'
            dotIndex = []
            #for row in range(len(tempData)):
            for index, char in enumerate(tempData[0]):
                if char == dot:
                    dotIndex.append(index)
            dotArray = []
            for row in range(len(tempData)):
                for dot in range(len(dotIndex)):
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'toneTime'] = np.array(altDotArray)         

        # !! ADD BEAM BREAKS !!
        itemIndex = 0
        found = False
        for line, row in enumerate(MPCArray):
            if 'M:' in row:
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
                    if tempData[row][dotIndex[dot]-5]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-4:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-4]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-3:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-3]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-2:dotIndex[dot]+4])
                    elif tempData[row][dotIndex[dot]-2]==' ':
                        dotArray.append(tempData[row][dotIndex[dot]-1:dotIndex[dot]+4])
            altDotArray = [float(item) for item in dotArray]
            dataSheet.at[0,'beamBreakTime'] = np.array(altDotArray)        
            
        dataSheet.to_excel('{}\\{}.xlsx'.format(dataPath,txtFiles[0][:-4]), index=True)
        return dataSheet