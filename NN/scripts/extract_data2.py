#!/home/amedina/build2/bin/python                                                                                                                                                   
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
from icecube.recclasses import I3LaputopParams, LaputopParameter
from icecube.stochastics import *
import multiprocessing as mp

## Are you using tableio?                                                                                                                                                                 
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

from scipy.optimize import curve_fit

import numpy as np
import pandas as pd
from glob import glob

directory = '/data/user/amedina/CosmicRay/Analysis_data/'
my_list = glob(directory + '*.i3.bz2')
output_file = '/data/user/amedina/CosmicRay/csv_files_data/data.csv'

def gh_function(depth,Xo,Xmax,lambda_value):
    x = (depth-Xo)/lambda_value
    m = (Xmax-Xo)/lambda_value
    n = ((x/m)**m) * np.exp(m-x)
    return np.log(n)

def get_chi2(x_true,x_predicted):
    x_predicted = np.array(x_predicted)
    x_true = np.array(x_true)
    chi2 = np.sum([((i-j)**2)/abs(j) for i,j in zip(x_true,x_predicted)])
    return chi2

def process_files(input_file):
    our_map = {}

    zenith = []
    S125 = []
    energy_loss = []
    he_stoch = []
    he_stoch2 = []
    
    A = []

    m_125 = []
    m_s2 = []
    m_s = []
    m_r = []
    m_o = []
    m_chi2 = []
    fit_status_m = []

    s_r = []
    s_o = []
    s_chi2 = []
    fit_status_s = []

    s_mean = []
    s_std = []

    charge = []
    N = []

    waveform_weight = []
    
    l3_file = dataio.I3File(input_file,'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics() #This gets the frame

        check = l3_fr['IT73AnalysisInIceQualityCuts'].values()
        if np.sum(check)/len(check) != 1:
            continue

        if 'Millipede_dEdX' in l3_fr:
            energy_loss.append(l3_fr['Stoch_Reco'].eLoss_1500)
            he_stoch.append(l3_fr['Stoch_Reco'].nHEstoch)
            he_stoch2.append(l3_fr['Stoch_Reco2'].nHEstoch)
        else:
            energy_loss.append(None)
            he_stoch.append(None)
            he_stoch2.append(None)

        waveform_weight.append(l3_fr['IceTopWaveformWeight'].value)

        zenith.append(l3_fr['Laputop'].dir.zenith)
        S125.append(l3_fr['LaputopParams'].value(LaputopParameter.Log10_S125))
        m_125.append(l3_fr['m_fit']['m_125'])
        m_o.append(l3_fr['m_fit']['m_o'])
        m_r.append(l3_fr['m_fit']['m_r'])
        m_s.append(l3_fr['m_fit']['m_s'])
        m_s2.append(l3_fr['m_fit']['m_s2'])
        m_chi2.append(l3_fr['m_fit']['chi2'])
        fit_status_m.append(l3_fr['m_fit']['fit_status'])

        s_r.append(l3_fr['s_fit']['s_r'])
        s_o.append(l3_fr['s_fit']['s_o'])
        s_chi2.append(l3_fr['s_fit']['chi2'])
        fit_status_s.append(l3_fr['s_fit']['fit_status'])

        if 's_mean' in l3_fr['m_fit'].keys():
            s_mean.append(l3_fr['m_fit']['s_mean'])
            s_std.append(l3_fr['m_fit']['s_std'])
        else:
            s_mean.append(0)
            s_std.append(0)


        count = 0
        qtot = 0
        for i in l3_fr['LaputopHLCVEM'].keys():
            if np.isnan(l3_fr['LaputopHLCVEM'][i][0].charge):
                continue
            else:
                count+=1
                qtot+=l3_fr['LaputopHLCVEM'][i][0].charge
            
        charge.append(qtot)
        N.append(count)
        
        #my_variables, beta is not related to the age of the shower
        
        A.append(l3_fr['CurvatureOnlyParams'].value(LaputopParameter.CurvParabA))
        
    l3_file.close()
    our_map['zenith'] = zenith
    our_map['S125'] = S125
    our_map['energy_loss'] = energy_loss
    our_map['he_stoch'] = he_stoch
    our_map['he_stoch2'] = he_stoch2
    our_map['A'] = A
    our_map['m_125'] = m_125
    our_map['m_r'] = m_r
    our_map['m_s'] = m_s
    our_map['m_s2'] = m_s2
    our_map['m_o'] = m_o
    our_map['m_chi2'] = m_chi2
    our_map['fit_status_m'] = fit_status_m
    our_map['s_r'] = s_r
    our_map['s_o'] = s_o
    our_map['s_chi2'] = s_chi2
    our_map['s_mean'] = s_mean
    our_map['s_std'] = s_std
    our_map['fit_status_s'] = fit_status_s
    our_map['charge'] = charge
    our_map['N'] = N
    our_map['waveform_weight'] = waveform_weight


    return our_map

count = 0

for i in my_list:
    print(i)
    try:
        if count == 0:
            new_map = process_files(i)
        else:
            event = process_files(i)
            for id,value in event.items():
                new_map[id]+=value 
        count+=1
    except:
        continue

df = pd.DataFrame(new_map)
df.to_csv(output_file)
