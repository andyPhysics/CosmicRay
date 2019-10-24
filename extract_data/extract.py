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
from multiprocessing import Pool
from functools import partial
## Create ability to change settings from terminal ##                                                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",type=str,default='Level5_IC86.2013_genie_numu.014640.00000?.i3.bz2',
                    dest="input_file", help="name of the input file")
parser.add_argument("-n", "--name",type=str,default='Level5_IC86.2013_genie_numu.014640.00000X',
                    dest="output_name",help="name for output file (no path)")
parser.add_argument("-m","--mass",type=float,default=1.0,
                    dest="mass",help="This is the mass of the input file events")

args = parser.parse_args()
output_name = args.output_name

directory = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/'

my_files = get_file_list([directory+'12360/',directory+'12632/',directory+'12634/',directory+'12636/'])
#output = read_root_files(my_files,1)

p = Pool(5)
function = partial(read_root_files,input_mass=1)
output = p.map(function,my_files)

np.save(output_name,output)
