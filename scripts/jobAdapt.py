#!/usr/bin/python
#PBS -N Adapter
#PBS -V 
#PBS -l walltime=999:00:00,file=20000000kb

import sys
import os
import subprocess
import glob

#---------------------------------------------------------------
#Script to run through data files creating weight files
#-----------------------------------------------------------------

seed        = str(int(os.environ['Seed'])+int(os.environ['PBS_ARRAYID']))
print(seed)

outFolder   = os.environ['OutDir']
outFile     = "Events_"+str(seed)+".root"

if os.path.exists(outFolder + outFile):
    sys.exit("file exists")
    
outMain     = os.environ['OutMain']

#Static variables 
number      = os.environ['Number']   # Number of events
beamFile    = os.environ['beamSettings']
targetFile  = os.environ['radiatorSettings']

adapterBase = os.environ['sampler']

adapterFile = adapterBase+".astate"
adapterOut  = adapterBase+str(seed)+".astate"

#Directorys
ProcessingDir = "/scratch/"os.environ['USER']"/"
if not os.path.exists(ProcessingDir):
    os.makedirs(ProcessingDir)

if os.path.exists(outMain+adapterFile):
    subprocess.call(["cp", outMain+adapterFile, ProcessingDir+adapterFile])

subprocess.call(["/home/simong/Generators/GlueXCB/CBGen", "-d", ProcessingDir, "-f", outFile, "-s", seed, "-n", number, "-b", beamFile, "-t", targetFile, "-a", adapterFile])

# Copy output events and adapter stats to work disks 
subprocess.call(["cp", ProcessingDir+outFile, outFolder+outFile])
subprocess.call(["cp", ProcessingDir+adapterOut, outFolder+adapterOut])

# Remove working files
subprocess.call(["rm", ProcessingDir+outFile])
subprocess.call(["rm", ProcessingDir+adapterOut])
subprocess.call(["rm", ProcessingDir+adapterFile])

print "DONE"
