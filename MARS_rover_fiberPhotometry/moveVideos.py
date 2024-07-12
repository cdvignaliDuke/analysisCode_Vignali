# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:05:25 2024

@author: Carlo Vignali

@def: This will run through all video files and put them in their associated animal\\day folder

"""
#%% Part 0.0 - import packages
import os
import time
import subprocess
os.chdir('S:\\Private\\Data\\Vignali cv105\\code\\MARS')
import localFunctions as lf
os.environ["PATH"] += os.pathsep + 'C:\\Program Files\\ffmpeg\\bin\\'

#%% Part 1.0 - iterate through each mouse
inputPath = ('S:\\Private\\Data\\Vignali cv105\\data\\MARS\\')
storage = ('F:\\MARS\\videos\\')
allMice = [file for file in os.listdir(inputPath) if not os.path.isfile(os.path.join(inputPath, file))]

for m in range(len(allMice)):
    print('startTime for mouse {mus}: {t}'.format(mus=allMice[m],t=time.strftime("%H:%M:%S", time.localtime())))
    
#% Part 1.1 - compress .avi into .mp4 with HEVC codex    
    toCompress = [file for file in os.listdir(inputPath) if 'T0{}'.format(allMice[m]) in file and file.endswith('.avi') and os.path.isfile(os.path.join(inputPath, file))]
    
    if not toCompress == []:
        for vid in range(len(toCompress)):
            command = [
            'ffmpeg',
            '-i', os.path.join(inputPath, toCompress[vid]), # Input file path
            '-c:v', 'libx265', # Specify video codec
            '-crf', '17', # Constant Rate Factor (quality-level, where a lower number is higher quality, 17-28 good for HD)
            '-preset', 'slow', # Speed/quality tradeoff (veryfast, faster, fast, medium, slow)
            '-r', '200', # Output framerate
            os.path.join(inputPath, os.path.splitext(toCompress[vid])[0] + '.mp4') # Output file path
        ]
            try:
                subprocess.run(command, check=True)
                print("Compression complete.")
            except subprocess.CalledProcessError:
                print("Error during compression.")
          

#% Part 1.2 - move .avi into storage and .mp4 into data folder of mouse    
        lf.moveFiles(allMice[m],inputPath,prefix='T0',fileType='.mp4')
        lf.moveFiles(allMice[m],inputPath,outputPath=storage,prefix='T0',fileType='.avi')
    
        print('endTime for mouse {mus}: {t}'.format(mus=allMice[m],t=time.strftime("%H:%M:%S", time.localtime())))

