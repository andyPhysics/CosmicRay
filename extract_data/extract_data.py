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
from scipy.optimize import curve_fit
from scipy.stats import chisquare


## Create ability to change settings from terminal ##

def Gaisser_hillas_function(x,m,alpha,b,a):
    n = m*np.log(alpha*(x-a))-alpha*x + b
    return n

def Gaisser_exp(x,m,alpha,b,a):
    n = np.exp(Gaisser_hillas_function(x,m,alpha,b,a))
    return n

def get_file_list(directories):
    file_list = []
    for i in directories:
        files = os.listdir(i)
        for j in files:
            file_list.append(i + j)
    return file_list

def get_Xmax(depth,num):
    x = zip(depth,num)
    output = []
    for i in x:
        if i[1] > np.exp(-8):
            output.append(i)
        else:
            continue
    length = len(output)
    depth_new = np.array(list(zip(*output))[0][0:length-3])
    num_new = np.array(list(zip(*output))[1][0:length-3])
    popt,pcov = curve_fit(Gaisser_exp,depth_new,num_new,bounds=((0,0,-np.inf,-np.inf),(np.inf,np.inf,np.inf,min(depth))))
    return popt,depth_new,num_new

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
    run = []
    event = []
    mass = []
    energy = []
    s70 = []
    s150 = []
    s125 = []
    beta = []
    zenith = []
    azimuth = []
    x1 = []
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
    stoch_energy2 = []
    rel_stoch_energy2 = []
    chi2_red2 = []
    stoch_depth2 = []
    n_he_stoch2 = []
    fit_status2 = []
    mc_weight = []
    num_EPlus = []
    num_EMinus = []
    num_ETotal = []
    depth = []
    A = []
    D = []
    N = []
    chi2_curvature = []

    for i in files:
        x = uproot.open(i)
        if len(x.keys()) == 0 :
            continue
        Stoch_Reco_Found = False
        for f in x.keys():
            if f=='Stoch_Reco;1':
                Stoch_Reco_Found = True
        if not Stoch_Reco_Found:
            continue
        
        depth += [x['MCPrimaryInfo']['longDepth'].array()]
        num_EPlus += [x['MCPrimaryInfo']['longNumEMinus'].array()]
        num_EMinus += [x['MCPrimaryInfo']['longNumEPlus'].array()]
        run += [x['I3EventHeader']['Run'].array()]
        event += [x['I3EventHeader']['Event'].array()]
        mass += [input_mass]*len(depth)
        energy += [x['MCPrimary']['energy'].array()]
        s70 += [x['LaputopParams']['s70'].array()]
        s150 += [x['LaputopParams']['s150'].array()]
        s125 += [x['LaputopParams']['s125'].array()]
        beta += [x['LaputopParams']['beta'].array()]
        zenith += [x['MCPrimary']['zenith'].array()]
        azimuth += [x['MCPrimary']['azimuth'].array()]
        x1 += [x['MCPrimary']['x'].array()]
        y += [x['MCPrimary']['y'].array()]
        z += [x['MCPrimary']['z'].array()]
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
        stoch_energy2 += [x['Stoch_Reco2']['stoch_energy'].array()]
        rel_stoch_energy2 += [x['Stoch_Reco2']['rel_stoch_energy'].array()]
        chi2_red2 += [x['Stoch_Reco2']['chi2_red'].array()]
        stoch_depth2 += [x['Stoch_Reco2']['stoch_depth'].array()]
        n_he_stoch2 += [x['Stoch_Reco2']['n_he_stoch'].array()]
        fit_status2 += [x['Stoch_Reco2']['fit_status'].array()]
        mc_weight += [x['MCPrimaryInfo']['weight'].array()]
        A += [x['CurvatureOnlyParams']['A'].array()]
        D += [x['CurvatureOnlyParams']['D'].array()]
        N += [x['CurvatureOnlyParams']['N'].array()]
        chi2_curavture += [x['CurvatureOnlyParams']['chi2'].array()]
        
    my_dict = dict(run = run,
                   event = event,
                   mass = mass,
                   num_EPlus = num_EPlus,
                   num_EMinus = num_EMinus,
                   depth = depth,
                   energy = energy,
                   s70 = s70,
                   s150 = s150,
                   s125 = s125,
                   beta = beta,
                   zenith = zenith,
                   azimuth = azimuth,
                   x = x1,
                   y = y,
                   z = z,
                   eloss_1500 = eloss_1500,
                   eloss_1800 = eloss_1800,
                   eloss_2100 = eloss_2100,
                   eloss_2400 = eloss_2400,
                   stoch_energy = stoch_energy,
                   rel_stoch_energy = rel_stoch_energy,
                   chi2 = chi2,
                   chi2_red = chi2_red,
                   stoch_depth = stoch_depth,
                   n_he_stoch = n_he_stoch,
                   fit_status = fit_status,
                   stoch_energy2 = stoch_energy2,
                   rel_stoch_energy2 = rel_stoch_energy2,
                   chi2_red2 = chi2_red2,
                   stoch_depth2 = stoch_depth2,
                   n_he_stoch2 = n_he_stoch2,
                   fit_status2 = fit_status2,
                   mc_weight = mc_weight,
                   A = A,
                   D = D,
                   N = N,
                   chi2_curvature = chi2_curvature)
    return my_dict

