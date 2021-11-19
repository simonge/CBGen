#!/usr/bin/python3

import os
import sys
import subprocess
import glob
import time

#-----------------------------------------------------------------
#Script to run through data files creating weight files
#-----------------------------------------------------------------

    


#Name of the script to run
scriptNames  = ["/home/simong/Generators/GlueXCB/scripts/jobCBGen.py","/home/simong/Generators/GlueXCB/scripts/jobAdapter.py"]
jobBaseNames = ["CBGen","Adapter"]
nFiles       = [200,1] #100
#nFiles       = [0,1] #100
#nFiles       = [1000,1] #100

outBase = "/w/work5/home/simong/CBremGen/Adapt/"
os.environ['OutBase'] = outBase

seedStart        = 21101
nStart           = 1
nRepeats         = 10 #400
#nStart           = 104
#nRepeats         = 111 #400
outputROOT       = [0,10,20,30,40,50,60,62,70,80,100,110]
maxq             = 200
number           = 1000000
#number           = 500000
#beamSettings     = "/home/simong/Generators/GlueXCB/setup/MainzSettingsFara.txt"
beamSettings     = "/home/simong/Generators/GlueXCB/setup/MainzSettingsNew.txt"
#beamSettings     = "/home/simong/Generators/GlueXCB/setup/MainzSettings.txt"
#beamSettings     = "/home/simong/Generators/GlueXCB/setup/GlueXSettings.txt"
#radiatorSettings = "/home/simong/Generators/GlueXCB/setup/Diamond.txt"
#radiatorSettings = "/home/simong/Generators/GlueXCB/setup/DiamondNew.txt"
base             = "Diamond60mrad"
#radiatorSettings = "/home/simong/Generators/GlueXCB/setup/Copper.txt"
#base             = "Copper3"
adaptBase        = "sampler"

    
outDirectory = outBase + base + "/"
if not os.path.exists(outDirectory):
    os.makedirs(outDirectory)
    
os.environ['OutMain'] = outDirectory    

os.environ['sampler'] = adaptBase
    
os.environ['radiatorSettings'] = radiatorSettings
os.environ['beamSettings']     = beamSettings
os.environ['Number']           = str(number)
os.environ['OutputROOT']       = str(0)
    
jobprocessed = 0
jobcount     = 0
previousJob  = "none"

seedStart += (nStart-1)*sum(nFiles)

#Loop over adaptations
for repeat in range(nStart,nRepeats):

    os.environ['repeat']  = str(repeat)
    
    repeatDir = outDirectory + "Repeat_" + str(repeat) + "/"

    if not os.path.exists(repeatDir):
        os.makedirs(repeatDir)

    os.environ['OutDir']  = repeatDir
    os.environ['OutputROOT'] = str(int(repeat in outputROOT))
    print(repeat,os.environ['OutputROOT'])
    
    # Add job array to the queue instead of a loop    
    for script, job, files, in zip(scriptNames,jobBaseNames,nFiles):
        os.environ['Seed'] = str(seedStart)
        # Range of job numbers
        JobRange = "1-"+str(files)+"%"+str(maxq)
        thisJob  = job + str(repeat)
        
        # Sumbit jobs
        if(previousJob=="none"):
            process = subprocess.run(["qsub","-q","jarray", "-N", thisJob, "-t", JobRange, script],env=dict(os.environ),stdout=subprocess.PIPE,universal_newlines=True)
        else:
            process = subprocess.run(["qsub","-q","jarray", "-N", thisJob, "-t", JobRange, "-W", "depend=afteranyarray:"+previousJob, script],env=dict(os.environ),stdout=subprocess.PIPE,universal_newlines=True)
        previousJob = process.stdout.rstrip()
        print(previousJob)
        print("HI")
        seedStart  += files
    

print('All Done')


