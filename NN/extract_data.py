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


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",type=str,default='12360',
                    dest="directory_number",help="directory number")
args = parser.parse_args()

def process_files(file_tuple):
    print(file_tuple)
    outdir = '/data/user/amedina/CosmicRay/I3_Files/%s'%(args.directory_number)

    assert os.path.isfile(file_tuple[1]), "file doesn't exist: {}".format(file_tuple[1])
    assert os.path.isfile(file_tuple[0]), "file doesn't exist: {}".format(file_tuple[0])
    l3_file = dataio.I3File(file_tuple[1],'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics() #This gets the frame

    l3_file.close()


