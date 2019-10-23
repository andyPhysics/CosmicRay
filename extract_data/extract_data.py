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
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.optimize import brentq

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
    popt,pcov = curve_fit(Gaisser_hillas_function,depth,num,bounds=((0,0,-np.inf,-np.inf),(np.inf,np.inf,np.inf,min(depth))),p0=[1,1,-1,10])
    return popt

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
    lambda_values = []
    X_o = []
    chi2_xmax = []
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
    time_start_mjd_sec = []
    time_start_mjd_ns = []
    time_start_mjd_day = []
    stoch_energy2 = []
    rel_stoch_energy2 = []
    chi2_red2 = []
    eloss_1500_red = []
    stoch_depth2 = []
    n_he_stoch2 = []
    fit_status2 = []
    mc_weight = []
    sum_value_prediction = []
    depth_reduced = []
    depth2 = []
    sum_value2 = []
    
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
        sum_value = [i/max(i) for i in sum_value]
        depth2.append(depth)
        sum_value2.append(sum_value)
        count = 0
        for i in range(depth.shape[0]):
            count+=1
            new_values = zip(depth[i],sum_value[i])
            new_values2 = []
            for j in new_values:
                if j[1] < np.exp(-6):
                    continue
                else:
                    new_values2.append(j)
            depth1 = np.array(list(zip(*new_values2))[0])
            sum_value1 = np.log(list(zip(*new_values2))[1])
            prediction = get_Xmax(depth1,sum_value1)
            xmax.append(prediction[0]/prediction[1]+prediction[3])
            lambda_values.append(1/prediction[1])
            X_o.append(prediction[3])
  #          X_o.append(brentq(Gaisser_exp,1e-100,prediction[0]/prediction[1]+prediction[3],args=(prediction[0],prediction[1],prediction[2],prediction[3])))
            chi2_xmax.append(chisquare(Gaisser_exp(depth1,prediction[0],prediction[1],prediction[2],prediction[3]),f_exp=list(zip(*new_values2))[1],ddof=4)[0])
            sum_value_prediction.append(Gaisser_exp(depth1,prediction[0],prediction[1],prediction[2],prediction[3]))
            depth_reduced.append(depth1)


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
        mc_weight += [x['MCPrimaryInfo']['weight'].array()]


    my_dict = dict(run = np.hstack(run),
                   event = np.hstack(event),
                   mass = np.hstack(mass),
                   energy = np.hstack(energy),
                   depth = np.hstack(depth2),
                   sum_value = np.hstack(sum_value2),
                   xmax = np.hstack(xmax),
                   lambda_values = np.hstack(lambda_values),
                   X_o = np.hstack(X_o),
                   chi2_xmax = np.hstack(chi2_xmax),
                   sum_value_prediction = np.array(sum_value_prediction),
                   depth_reduced = np.array(depth_reduced),
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

