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

def mergeDict(dict1, dict2):
   dict3 = {}
   keys = dict1.keys()
   for i in keys:
      print(dict1[i])
      dict3[i] = np.hstack([dict1[i],dict2[i]])
   return dict3

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name",type=str,default='Test.npy',
                    dest="output_name",help="name for output file (no path)")
parser.add_argument("-m","--mass",type=float,default=1.0,
                    dest="mass",help="This is the mass of the input file events")
parser.add_argument("-f","--files",type=str,nargs='+',
                     dest="files",help="This is the list of directories as numbers such as 123456/")

args = parser.parse_args()
output_name = args.output_name

directory = '/data/user/amedina/CosmicRay/Curvature/'
output_directories = []
print(args.files)
for i in args.files:
    output_directories.append(directory+i+'/')

my_files = get_root_files(output_directories)

shared_dict,InIce_cuts,IceTop_cuts = read_root_files(my_files,args.mass)

First_dict = InIce_cuts[0]
for i in InIce_cuts[1:]:
    First_dict = mergeDict(First_dict,i)

Second_dict = IceTop_cuts[0]
for i in IceTop_cuts[1:]:
    Second_dict = mergeDict(Second_dict,i)

print(Second_dict)

shared_dict['InIce_cuts']=First_dict
shared_dict['IceTop_cuts']=Second_dict
print(shared_dict.keys())
np.save(output_name,shared_dict)
