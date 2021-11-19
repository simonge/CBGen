#!/usr/bin/python

import os
import sys
import subprocess
import glob
import time

#---------------------------------------------------------------
#Script to run through data files creating weight files
#-----------------------------------------------------------------

outDirectory = "/w/work1/home/simong/CBremGen/Events/"
os.environ['OutDir'] = outDirectory

seedStart = 0
nFiles    = 100
maxq      = 200
number    = 1000000
beamE     = 1.6
cohEdge   = 0.45

os.environ['beamE'] = str(beamE)
os.environ['cohEdge'] = str(cohEdge)
os.environ['Number']   = str(number)
    
jobprocessed = 0

for i in range(0,nFiles):
    seed = seedStart+i
    os.environ['Seed'] = str(seed)
    
    outFile = "Events_"+str(seed)+".root"
    os.environ['OutFile'] = outFile
        
    if os.path.exists(outDirectory + outFile):
        continue
    
    jobprocessed+=1
    subprocess.call(["qsub", "/home/simong/Generators/GlueXCB/scripts/jobSimNew.py"],env=dict(os.environ))
    time.sleep(0.01)
    jobcount = subprocess.check_output(["qstat | grep \'R\|Q\'"], shell=True).count(os.environ['USER'])
    print jobcount
    while jobcount>maxq:
        time.sleep(1)
        jobcount = subprocess.check_output(["qstat | grep \'R\|Q\'"], shell=True).count(os.environ['USER'])



print 'All Done'

