#!/usr/bin/python
#PBS -N StoneHengeTest
#PBS -V 
#PBS -l walltime=999:00:00,file=20000000kb

import os
import subprocess
import glob

#---------------------------------------------------------------
# Script to run through data files creating weight files
#-----------------------------------------------------------------

seed      = '52'
outFolder = '/scratch/simong/'
outFile   = 'out.root'
analysis  = '2'
number    = '350'
colR      = '0.0025'
colRMS    = '0.00034'
colD      = '3'
beamEmm   = '25e-9'
thickness = '300e-6'
thickCu   = '100e-6'
yOff      = '5.00'
beamE     = '1.508'

subprocess.call(['/home/simong/Generators/CohBremRJ/build/CohGen', '-f', outFolder, '-s', seed, '-a', analysis, '-n', number, '-c', colR, '-d', colD, '-k', colRMS, '-b', beamEmm, '-t', thickness, '-y', yOff, '-u', thickCu, '-e', beamE])

print 'DONE'
