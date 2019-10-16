#!/usr/bin/env python

import numpy as np
import h5py
import argparse

from icecube import icetray, dataio, dataclasses, simclasses, recclasses
from icecube.recclasses import LaputopParameter as Par
from I3Tray import I3Units

from collections import OrderedDict
import itertools
import random
import datetime
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

def get_Xmax(depth,num):
    max_num = max(num)
    index = 0
    for i in num:
        if i == max_num:
            break
        else:
            i+=1
    depth_cut = [depth[index-1],depth[index],depth[index+1]]
    num_cut = [num[index-1],num[index],num[index+1]]
    x = np.polyfit(depth_cut,num_cut,deg=2)
    max_value = -x[1]/(2*x[0])
    return max_value

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
        xmax = []
        run = []
        event = []
        s125 = []
        beta = []
        zenith = []
        azimuth = []
        x = []
        y = []
        z = []
        eloss_1500 = []
        eloss_1800 = []
        eloss_2100 = []
        eloss_2400 = []
        stoch_energy = []
        rel_stoch_energy = []
        chi2 = []
        chi2_red = []
        stoch_depth = []
        n_he_stoch = []
        fit_status = []
        time_start_mjd_sec = []
        time_start_mjd_ns = []
        time_start_mjd_day = []
        eloss_1500_red = []
        stoch_energy2 = []
        rel_stoch_energy2 = []
        chi2_red2 = []
        stoch_depth2 = []
        n_he_stoch2 = []
        fit_status2 = []
        eloss_1500_red2 = []
        mc_weight = []
        nch = []
        qtot = []
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
            sum_values = np.array(event_numEMinus)+np.array(event_numEPlus)
            Xmax = xmax.append(get_Xmax(event_depth,sum_values))
            depth.append(event_depth)
            mass.append(args.mass)
            energy.append(frame['MCPrimary'].energy)
            run.append(frame['I3EventHeader'].run_id)
            event.append(frame['I3EventHeader'].event_id)
            s125.append(frame['LaputopSmallParams'].value(Par.Log10_S125))
            beta.append(frame['LaputopSmallParams'].value(Par.Beta))
            zenith.append(frame['MCPrimary'].dir.zenith)
            azimuth.append(frame['MCPrimary'].dir.azimuth)
            x.append(frame['MCPrimary'].pos.x)
            y.append(frame['MCPrimary'].pos.y)
            z.append(frame['MCPrimary'].pos.z)
            fit_status.append(frame['LaputopSmall'].fit_status)
            time_start_mjd_sec.append(frame['I3EventHeader'].start_time)
            print(str(frame['I3EventHeader'].start_time).split()[1].replace(',',''))

#            print(frame['IsSmallShower'].value)
        print(xmax)
        my_dict = dict(numEMinus=numEMinus,
                       numEPlus=numEPlus,
                       depth=depth,
                       mass=mass,
                       energy=energy,
                       xmax=xmax)
        # close the input file once we are done
        del event_file
    return my_dict
my_files = ['/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12360/Level3_IC86.2012_12360_Run006667.i3.gz']


values = read_files(my_files)

np.save('First.npy',values)
