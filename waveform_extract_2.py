#!/home/amedina/build_stable/bin/python                                                                                                                                                   
import argparse
import os,sys,getopt
#import logging                                                                                                                                                                           
from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, frame_object_diff,simclasses,WaveCalibrator,tpx
from icecube import gulliver, gulliver_modules, lilliput,photonics_service,wavedeform,wavereform
from icecube.payload_parsing import I3DOMLaunchExtractor
load('millipede')
load('stochastics')

from icecube.icetray import I3Module
from icecube.dataclasses import I3EventHeader, I3Particle
from icecube.recclasses import I3LaputopParams
from icecube.stochastics import *
import multiprocessing as mp

## Are you using tableio?                                                                                                                                                                 
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

import numpy as np


#level2 files are separated into 1000 per sub directory
level2_directory = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/12360/level2/0000000-0000999/' 
level2_file = level2_directory + 'Level2_IC86_corsika_icetop.010410.000002.i3.bz2'

#All files are in one directory for level3 files
level3_directory = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12360/'
level3_file = level3_directory + 'Level3_IC86.2012_12360_Run000002.i3.gz'


class Check_Laputop(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Configure(self):
        pass

    def Physics(self, frame):
        self.PushFrame(frame)

tray = I3Tray()

tray.AddModule("I3Reader","reader")(
        ("FileNameList", [level3_file])
    )

tray.Add("Dump")

tray.Execute(10)

