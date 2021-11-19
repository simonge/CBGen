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
jobNo = int(os.environ['PBS_ARRAYID'])
seed        = str(int(os.environ['Seed'])+jobNo)
repeat      = os.environ['repeat']
print(seed)

outFolder   = os.environ['OutDir']
outFile     = "Events_"+repeat+"_"+str(jobNo)+".root"

if os.path.exists(outFolder + outFile):
    sys.exit("file exists")
    
outMain     = os.environ['OutMain']
outputROOT  = os.environ['OutputROOT']
#Static variables 
number      = os.environ['Number']   # Number of events
beamFile    = os.environ['beamSettings']
targetFile  = os.environ['radiatorSettings']

adapterBase = os.environ['sampler']

if repeat==0:
    adapterFile = adapterBase+".astate"
else:
    adapterFile = adapterBase+"_"+str(int(repeat)-1)+".astate"
    
adapterOut  = adapterBase+"_"+repeat+"_"+str(jobNo)+".astate"

if os.path.exists(outFolder + adapterOut):
    sys.exit("adapter exists")

#Directorys
ProcessingDir = "/scratch/"+os.environ['USER']+"/temp"+seed+"/"
if not os.path.exists(ProcessingDir):
    os.makedirs(ProcessingDir)

if os.path.exists(outMain+adapterFile):
    subprocess.call(["cp", outMain+adapterFile, ProcessingDir+adapterFile])

subprocess.call(["/home/simong/Generators/GlueXCB/CBGen", "-d", ProcessingDir, "-f", outFile, "-s", seed, "-n", number, "-b", beamFile, "-t", targetFile, "-a", adapterFile, "-o", adapterOut, "-root", outputROOT ])

#subprocess.call(["ls", ProcessingDir,"-l"])
# Copy output events and adapter stats to work disks
if(outputROOT):
    subprocess.call(["cp", ProcessingDir+outFile, outFolder+outFile])
    #subprocess.call(["rm", ProcessingDir+outFile])

    
subprocess.call(["cp", ProcessingDir+adapterOut, outFolder+adapterOut])

# Remove working files
subprocess.call(["rm", '-r', ProcessingDir])
#subprocess.call(["rm", ProcessingDir+adapterFile])

print "DONE"
