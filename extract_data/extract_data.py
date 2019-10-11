#!/usr/bin/env python

import numpy as np
import h5py
import argparse

from icecube import icetray, dataio, dataclasses, simclasses
from I3Tray import I3Units

from collections import OrderedDict
import itertools
import random

import sys,os

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


def read_files(filename_list):

    output_features_IC = []
    output_labels = []

    for event_file_name in filename_list:
        print("reading file: {}".format(event_file_name))
        event_file = dataio.I3File(event_file_name)
        
        event_number = 0

        numEMinus = []
        numEPlus = []
        depth = []
        mass = []
        energy = []

        while event_file.more():
            frame = event_file.pop_physics()
            long_profile = frame['MCPrimaryInfo'].longProfile
            event_numEMinus = []
            event_numEPlus = []
            event_depth = []
            for i in long_profile:
                event_numEMinus.append(i.numEMinus)
                event_numEPlus.append(i.numEPlus)
                event_depth.append(i.depth)
            numEMinus.append(event_numEMinus)
            numEPlus.append(event_numEPlus)
            depth.append(event_depth)
            mass.append(args.mass)
            energy.append(frame['MCPrimary'].energy)
        my_dict = dict(numEMinus=numEMinus,
                       numEPlus=numEPlus,
                       depth=depth,
                       mass=mass,
                       energy=energy)
        # close the input file once we are done
        del event_file
    return my_dict
my_files = ['/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12360/Level3_IC86.2012_12360_Run006667.i3.gz']


values = read_files(my_files)

np.save('First.npy',values)
