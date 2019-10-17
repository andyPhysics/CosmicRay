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

def get_file_list(directories):
    file_list = []
    for i in directories:
        files = os.listdir(i)
        for j in files:
            file_list.append(i + j)
    return file_list

def get_Xmax(depth,num):
    max_num = max(num)
    index = 0
    depth_cut = []
    num_cut = []
    for i,j in zip(depth,num):
        if j == max_num:
            depth_cut.append(i)
            num_cut.append(j)
            break
        index += 1
    if index == depth.shape[0]-1:
        max_value = depth[-1]
    else:
        num_cut.append(num[index-1])
        num_cut.append(num[index+1])
        depth_cut.append(depth[index-1])
        depth_cut.append(depth[index+1])
        x = np.polyfit(depth_cut,num_cut,deg=2)
        if x[0] == 0:
            print(max_num,num_cut)
        max_value = -x[1]/(2*x[0])
    return max_value

def read_xmax_from_i3_file(event_file_name):
    print("reading file: {}".format(event_file_name))
    event_file = dataio.I3File(event_file_name)
    numEMinus = []
    numEPlus = []
    run = []
    event = []
    xmax = []

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
        xmax.append(get_Xmax(event_depth,sum_values))
        run.append(frame['I3EventHeader'].run_id)
        event.append(frame['I3EventHeader'].event_id)
    return zip(xmax,run,event)

def read_root_files(files,input_mass):
    count = 0
    run = []
    event = []
    mass = []
    energy = []
    xmax = []
    s70 = []
    s150 = []
    s125 = []
    beta = []
    zenith = []
    azimuth = []
    x1 = []
    y = []
    z = []
    itsiz = []
    iisiz = []#??? Not sure what this is                                                                                                                                            
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
    #    eloss_1500_red2 = [] Can't find                                                                                                                                            
    mc_weight = []
    #    nch = [] Can't find                                                                                                                                                        
    #    qtot = [] Can't find  

    for i in files:
        x = uproot.open(i)
        run += [x['I3EventHeader']['Run'].array()]
        event += [x['I3EventHeader']['Event'].array()]
        mass += [input_mass]*run[count].shape[0]
        energy += [x['MCPrimary']['energy'].array()]
        depth = x['MCPrimaryInfo']['longDepth'].array()
        numEPlus = x['MCPrimaryInfo']['longNumEMinus'].array()
        numEMinus = x['MCPrimaryInfo']['longNumEPlus'].array()
        sum_value = np.array(numEPlus)+np.array(numEMinus)
        xmax += [get_Xmax(i,j) for i,j in zip(depth,sum_value)]
        s70 += [x['LaputopParams']['s70'].array()]
        s150 += [x['LaputopParams']['s150'].array()]
        s125 += [x['LaputopParams']['s125'].array()]
        beta += [x['LaputopParams']['beta'].array()]
        zenith += [x['MCPrimary']['zenith'].array()]
        azimuth += [x['MCPrimary']['azimuth'].array()]
        x1 += [x['MCPrimary']['x'].array()]
        y += [x['MCPrimary']['y'].array()]
        z += [x['MCPrimary']['z'].array()]
        itsiz = []#?? Not sure what this is
        iisiz = []#??? Not sure what this is
        eloss_1500 += [x['Stoch_Reco']['eloss_1500'].array()]
        eloss_1800 += [x['Stoch_Reco']['eloss_1800'].array()]
        eloss_2100 += [x['Stoch_Reco']['eloss_2100'].array()]
        eloss_2400 += [x['Stoch_Reco']['eloss_2400'].array()]
        stoch_energy += [x['Stoch_Reco']['stoch_energy'].array()]
        rel_stoch_energy += [x['Stoch_Reco']['rel_stoch_energy'].array()]
        chi2 += [x['Stoch_Reco']['chi2'].array()]
        chi2_red += [x['Stoch_Reco']['chi2_red'].array()]
        stoch_depth += [x['Stoch_Reco']['stoch_depth'].array()]
        n_he_stoch += [x['Stoch_Reco']['n_he_stoch'].array()]
        fit_status += [x['Stoch_Reco']['fit_status'].array()]
        time_start_mjd_sec = []
        time_start_mjd_ns = []
        time_start_mjd_day = []
        eloss_1500_red = []
        stoch_energy2 += [x['Stoch_Reco2']['stoch_energy'].array()]
        rel_stoch_energy2 += [x['Stoch_Reco2']['rel_stoch_energy'].array()]
        chi2_red2 += [x['Stoch_Reco2']['chi2_red'].array()]
        stoch_depth2 += [x['Stoch_Reco2']['stoch_depth'].array()]
        n_he_stoch2 += [x['Stoch_Reco2']['n_he_stoch'].array()]
        fit_status2 += [x['Stoch_Reco2']['fit_status'].array()]
        #    eloss_1500_red2 = [] Can't find
        mc_weight += [x['MCPrimaryInfo']['weight'].array()]
        #    nch = [] Can't find
        #    qtot = [] Can't find
        count += 1
                
    my_dict = dict(run = np.hstack(run),
                   event = np.hstack(event),
                   mass = np.hstack(mass),
                   energy = np.hstack(energy),
                   xmax = np.hstack(xmax),
                   s70 = np.hstack(s70),
                   s150 = np.hstack(s150),
                   s125 = np.hstack(s125),
                   beta = np.hstack(beta),
                   zenith = np.hstack(zenith),
                   azimuth = np.hstack(azimuth),
                   x = np.hstack(x1),
                   y = np.hstack(y),
                   z = np.hstack(z),
                   eloss_1500 = np.hstack(eloss_1500),
                   eloss_1800 = np.hstack(eloss_1800),
                   eloss_2100 = np.hstack(eloss_2100),
                   eloss_2400 = np.hstack(eloss_2400),
                   stoch_energy = np.hstack(stoch_energy),
                   rel_stoch_energy = np.hstack(rel_stoch_energy),
                   chi2 = np.hstack(chi2),
                   chi2_red = np.hstack(chi2_red),
                   stoch_depth = np.hstack(stoch_depth),
                   n_he_stoch = np.hstack(n_he_stoch),
                   fit_status = np.hstack(fit_status),
                   stoch_energy2 = np.hstack(stoch_energy2),
                   rel_stoch_energy2 = np.hstack(rel_stoch_energy2),
                   chi2_red2 = np.hstack(chi2_red2),
                   stoch_depth2 = np.hstack(stoch_depth2),
                   n_he_stoch2 = np.hstack(n_he_stoch2),
                   fit_status2 = np.hstack(fit_status2),
                   mc_weight = np.hstack(mc_weight))
    return my_dict

my_files = ['/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/12360/Level3_IC86.2012_12360_Part099.root']#,'/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/12360/Level3_IC86.2012_12360_Part098.root']
#my_files = get_file_list(['/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/rootfiles/12360/'])
output = read_root_files(my_files,1)

np.save('First.npy',output)
