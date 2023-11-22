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
    our_map = {
        'zenith': [], 'S125': [], 'energy_loss': [], 'he_stoch': [], 'he_stoch2': [],
        'A': [], 'm_125': [], 'm_r': [], 'm_s': [], 'm_s2': [], 'm_o': [], 'm_chi2': [],
        'fit_status_m': [], 's_r': [], 's_o': [], 's_mean': [], 's_std': [], 's_chi2': [], 'fit_status_s': [],
        'charge': [], 'N': [], 'waveform_weight': []
    }
    
    l3_file = dataio.I3File(input_file,'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics() #This gets the frame

        check = l3_fr['IT73AnalysisInIceQualityCuts'].values()
        if np.sum(check)/len(check) != 1:
            continue

        if 'Millipede_dEdX' in l3_fr:
            our_map['energy_loss'].append(l3_fr['Stoch_Reco'].eLoss_1500)
            our_map['he_stoch'].append(l3_fr['Stoch_Reco'].nHEstoch)
            our_map['he_stoch2'].append(l3_fr['Stoch_Reco2'].nHEstoch)
        else:
            our_map['energy_loss'].append(None)
            our_map['he_stoch'].append(None)
            our_map['he_stoch2'].append(None)

        our_map['A'].append(l3_fr['CurvatureOnlyParams'].value(LaputopParameter.CurvParabA))
        our_map['zenith'].append(l3_fr['Laputop'].dir.zenith)
        our_map['S125'].append(l3_fr['LaputopParams'].value(LaputopParameter.Log10_S125))
        our_map['m_125'].append(l3_fr['m_fit']['m_125'])
        our_map['m_o'].append(l3_fr['m_fit']['m_o'])
        our_map['m_r'].append(l3_fr['m_fit']['m_r'])
        our_map['m_s'].append(l3_fr['m_fit']['m_s'])
        our_map['m_s2'].append(l3_fr['m_fit']['m_s2'])
        our_map['m_chi2'].append(l3_fr['m_fit']['chi2'])
        our_map['fit_status_m'].append(l3_fr['m_fit']['fit_status'])

        our_map['s_r'].append(l3_fr['s_fit']['s_r'])
        our_map['s_o'].append(l3_fr['s_fit']['s_o'])
        our_map['s_chi2'].append(l3_fr['s_fit']['chi2'])
        our_map['fit_status_s'].append(l3_fr['s_fit']['fit_status'])

        our_map['s_mean'].append(l3_fr['m_fit'].get('s_mean', 0))
        our_map['s_std'].append(l3_fr['m_fit'].get('s_std', 0))

        count = 0
        qtot = 0
        
        for i in l3_fr['LaputopHLCVEM'].values():
            if not np.isnan(i[0].charge):
                count += 1
                qtot += i[0].charge

        our_map['charge'].append(qtot)
        our_map['N'].append(count)

        try:
            our_map['waveform_weight'].append(l3_fr['IceTopWaveformWeight'].value)
        except KeyError:
            our_map['waveform_weight'].append(0)

    l3_file.close()

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
