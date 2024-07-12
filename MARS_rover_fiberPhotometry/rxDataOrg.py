# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:26:28 2024

@author: Carlo Vignali

@def: This is data organization code for fiber photometry data

@source: 
"""

def loadRawData(userInputs, experiment, seshNumber, thisMouse, thisDay, currFolders, rxSite):

    #%% Part 0.0 - import packages
    
    import os
    import numpy as np # fundamental package for scientific computing, handles arrays and maths
    from sklearn.metrics import auc
    import matplotlib.pyplot as plt  # standard Python plotting library
    import scipy.stats as stats
    import glob
    import pickle
    import tdt # tdt library
    import sys
    
    #
    
    os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
    import localFunctions as lf
    import rxPreproc 
    
    #%% Part 0.1 - set path and other core variables
    
    # !! SET PATH !!
    inputPath = ('S:\\Private\\Data\\Vignali cv105\\data\\MARS\\')
    outputPath = ('S:\\Private\\Data\\Vignali cv105\\analysis\\MARS\\')
    dataPath = inputPath + thisMouse + '\\' + thisDay + '\\'
    bxPath = outputPath + thisMouse + '\\bx\\'
    rxPath = outputPath + thisMouse + '\\' + thisDay + '\\fiber\\'
    os.makedirs(os.path.dirname(rxPath), exist_ok=True)
    
    # !! DIRECT TO FIBER DATA !!
    toLoad = [file for file in os.listdir(dataPath)] #and os.path.isfile(os.path.join(inputPath, file))]
    for file in range(len(toLoad)):
        if os.path.isdir(os.path.join(dataPath,toLoad[file])):
            rxDataFolder = toLoad[file] + '\\'
            
    #%% Part 1.0 - load data
    
    # !! LOAD DATATANK !!
    print('Importing rx data: ')
    tankPath = os.path.join(dataPath, rxDataFolder)
    pickles = [file for file in os.listdir(tankPath) if file.endswith('.pickle')]
    
    #%% Part 1.1 - import photometry stream
    if bool([file for file in os.listdir(tankPath) if file.endswith('.f32')]):
        print('   Fiber data found, pulling now: ')
        # !! PULL CHANNEL DATA !!
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
         
        fibDataSaved = {'data405_LH':data405_LH, 'data465_LH':data465_LH, 'data560_LH':data560_LH, 'data405_RH':data405_RH, 'data465_RH':data465_RH, 'data560_RH':data560_RH}
    
    else:    
        print('   No saved fiber data, loading from raw data tank: ')
        dataTank = tdt.read_block(tankPath, export='interlaced', outdir=tankPath, prefix='iso') #read length of data folder (tankPath)
        # note: the export and outdir flags allow the read_block function to pull out individual data streams for each sensor, and saves files in tankPath as F32 files (these are 32 bit floating point files - commonly used to represent sound waves)
        # ! see print(tdt.read_block.__doc__) for more details !
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
         
        fibDataSaved = {'data405_LH':data405_LH, 'data465_LH':data465_LH, 'data560_LH':data560_LH, 'data405_RH':data405_RH, 'data465_RH':data465_RH, 'data560_RH':data560_RH}

    #%% Part 1.2 - import epocs 
    file_path = 'epocSaved.pickle'    
    if file_path in pickles:
        if os.path.getsize(os.path.join(tankPath, file_path))>0:
            print('   Saved epoc files, pulling now: ')
            with open(os.path.join(tankPath, file_path), 'rb') as file:
                epocSaved = pickle.load(file)
        elif os.path.getsize(os.path.join(tankPath, file_path))==0:
            print('   Epoc data corrupted, loading from raw data tank: ')
            with open(tankPath + 'iso_interlaced_export.txt','r') as dataLines:
                fiberInfo = dataLines.readlines()
                            
                # !! FIND FS !!
                itemIndex = 0
                found = False
                for line, row in enumerate(fiberInfo):
                    if 'Freq:' in row:
                        itemIndex = line
                        found = True
                        break
                if found:
                    fs = float(fiberInfo[itemIndex][7:-1])
                else:
                    print('No "Freq" found in file. Program will continue, but make note that file may be improperly formatted.')    
                
            #%% Part 1.3 - import behavioral epocs FROM TDT DATATANK 
        
            epocData = tdt.read_block(tankPath, evtype=['epocs'])
                            
            epoc_ticker = epocData.epocs.Tick.onset
            epoc_onsetCues = epocData.epocs.sial.onset #house light off [+2.5s] tone plays
            epoc_onsetITI = epocData.epocs.sITI.onset
            epoc_beamBreak = epocData.epocs.seam.onset
        
            # used Epoc marker as a duplicate 
            epoc_startTrial = epocData.epocs.EpT_.onset
            epoc_startITI = epocData.epocs.EpI_.onset
            epoc_startBeam = epocData.epocs.EpC_.onset
         
            #%% Part 1.4 - save epoc information 
            
            epocSaved = {'epoc_ticker':epoc_ticker, 'epoc_onsetCues':epoc_onsetCues, 'epoc_onsetITI':epoc_onsetITI, 'epoc_beamBreak':epoc_beamBreak, 'epoc_startTrial':epoc_startTrial, 'epoc_startITI':epoc_startITI, 'epoc_startBeam':epoc_startBeam, 'fs':fs}
            file_path = 'epocSaved.pickle'    
            # Writing the dictionary to a file using pickle
            with open(os.path.join(tankPath, file_path), 'wb') as file:
                pickle.dump(epocSaved, file)     
    else:
        print('   No saved epoc data, loading from raw data tank: ')
        with open(tankPath + 'iso_interlaced_export.txt','r') as dataLines:
            fiberInfo = dataLines.readlines()
                        
            # !! FIND FS !!
            itemIndex = 0
            found = False
            for line, row in enumerate(fiberInfo):
                if 'Freq:' in row:
                    itemIndex = line
                    found = True
                    break
            if found:
                fs = float(fiberInfo[itemIndex][7:-1])
            else:
                print('No "Freq" found in file. Program will continue, but make note that file may be improperly formatted.')    
            
        #%% Part 1.3 - import behavioral epocs FROM TDT DATATANK 
    
        epocData = tdt.read_block(tankPath, evtype=['epocs'])
                        
        epoc_ticker = epocData.epocs.Tick.onset
        epoc_onsetCues = epocData.epocs.sial.onset #house light off [+2.5s] tone plays
        epoc_onsetITI = epocData.epocs.sITI.onset; 
        if epoc_onsetITI[0] > 1: epoc_onsetITI = np.insert(epoc_onsetITI,0,[0.0]); #session begins with ITI
        epoc_beamBreak = epocData.epocs.seam.onset
    
        # used Epoc marker as a duplicate 
        epoc_startTrial = epocData.epocs.EpT_.onset
        epoc_startITI = epocData.epocs.EpI_.onset
        epoc_startBeam = epocData.epocs.EpC_.onset
     
        #%% Part 1.4 - save epoc information 
        
        epocSaved = {'epoc_ticker':epoc_ticker, 'epoc_onsetCues':epoc_onsetCues, 'epoc_onsetITI':epoc_onsetITI, 'epoc_beamBreak':epoc_beamBreak, 'epoc_startTrial':epoc_startTrial, 'epoc_startITI':epoc_startITI, 'epoc_startBeam':epoc_startBeam, 'fs':fs}
        file_path = 'epocSaved.pickle'    
        # Writing the dictionary to a file using pickle
        with open(os.path.join(tankPath, file_path), 'wb') as file:
            pickle.dump(epocSaved, file)
        
    #%% Part 1.5 - import behavioral epocs   
    
    dataName_indivFP = ['data405_LH','data465_LH','data560_LH','data405_RH','data465_RH','data560_RH']
    colors_indivPlot = ['#0000FF','#008000','#FF0000','#ADD8E6','#90EE90','#FFC0CB']
    for f in range(len(dataName_indivFP)):
        if eval(dataName_indivFP[f]).any(): 
            plt.figure(figsize=(10, 4))
            plt.plot(eval(dataName_indivFP[f])[500:], color=colors_indivPlot[f])
            plt.title('Time Series for {}'.format(dataName_indivFP[f]))
            plt.xlabel('approx. time (ms)')
            plt.ylabel('raw fluorescence')
            plt.savefig('{directory}rawFluor_{chan}_{hemi}_{exp}{sesh}_{mus}.pdf'.format(directory=rxPath, chan=dataName_indivFP[f], hemi=rxSite, exp=experiment, sesh=seshNumber, mus=thisMouse), format='pdf', bbox_inches='tight', dpi=300)
            plt.show()   
    #%% Part 3.0 - save/output
    print(' Done.')
    return fibDataSaved, epocSaved
