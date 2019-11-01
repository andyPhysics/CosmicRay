#!/usr/bin/env python                                                                                                                                                                                       
import numpy as np
import h5py
import argparse
import uproot
from collections import OrderedDict
import itertools
import random
import datetime
import sys,os
from extract_data import *
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
## Create ability to change settings from terminal ##                                                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name",type=str,default='Test.npy',
                    dest="output_name",help="name for output file (no path)")
parser.add_argument("-m","--mass",type=float,default=1.0,
                    dest="mass",help="This is the mass of the input file events")
parser.add_argument("-f","--files",type=str,nargs='+',
                     dest="files",help="This is the list of directories as numbers such as 123456/")

args = parser.parse_args()
output_name = args.output_name

#directory = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/'
output_directories = []
print(args.files)
for i in args.files:
    output_directories.append(directory+i+'/')

my_files = get_file_list(output_directories)

shared_dict = read_root_files(my_files,args.mass)

np.save(output_name,shared_dict)
