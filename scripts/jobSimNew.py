#!/usr/bin/python
#PBS -N CBGen
#PBS -V 
#PBS -l walltime=999:00:00,file=20000000kb

import os
import subprocess
import glob

#---------------------------------------------------------------
#Script to run through data files creating weight files
#-----------------------------------------------------------------

seed      = os.environ['Seed']
outFolder = os.environ['OutDir']
outFile   = os.environ['OutFile']
number    = os.environ['Number']
cohEdge   = os.environ['cohEdge']
beamE     = os.environ['beamE']

#Directorys
ProcessingDir = "/scratch/simong/"
if not os.path.exists(ProcessingDir):
    os.makedirs(ProcessingDir)

subprocess.call(["/home/simong/Generators/GlueXCB/CBGen", "-d", ProcessingDir, "-f", outFile, "-s", seed, "-n", number, "-e", beamE, "-p", cohEdge])

subprocess.call(["cp", ProcessingDir+outFile, outFolder+outFile])

subprocess.call(["rm", ProcessingDir+outFile])

print "DONE"
