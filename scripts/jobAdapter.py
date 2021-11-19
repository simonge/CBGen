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


inFolder     = os.environ['OutDir']
outFile      = os.environ['OutMain']+os.environ['sampler']+"_"+os.environ['repeat']+".astate"

if os.path.exists(outFile):
    sys.exit("sampler exists")
    
pro = subprocess.call(["/home/simong/Generators/GlueXCB/build/Adapt","-d",inFolder,"-o",outFile])


print "DONE"
