#!/usr/bin/env python                                                                                                                                                                                       

import numpy as np
import h5py
import argparse

from icecube import icetray, dataio, dataclasses, simclasses, recclasses
from icecube.recclasses import LaputopParameter as Par
from I3Tray import I3Units
import uproot
from collections import OrderedDict
import itertools
import random
import datetime
import sys,os
from extract_data import *
## Create ability to change settings from terminal ##                                                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",type=str,default='Level5_IC86.2013_genie_numu.014640.00000?.i3.bz2',
                    dest="input_file", help="name of the input file")
parser.add_argument("-n", "--name",type=str,default='Level5_IC86.2013_genie_numu.014640.00000X',
                    dest="output_name",help="name for output file (no path)")
parser.add_argument("-m","--mass",type=float,default=1.0,
                    dest="mass",help="This is the mass of the input file events")

args = parser.parse_args()
input_file = args.input_file
output_name = args.output_name

my_files = ['/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/12360/Level3_IC86.2012_12360_Part099.root']

output = read_root_files(my_files,1)

np.save('First.npy',output)
