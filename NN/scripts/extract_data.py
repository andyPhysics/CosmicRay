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

directory = '/data/user/amedina/CosmicRay/Analysis/'
number = sys.argv[1]
my_list = glob(directory + '/%s_*.i3.bz2'%(number))
output_file = sys.argv[2]
print(my_list)

def gh_function(depth,Xo,Xmax,lambda_value):
    #This is the Gaisser-Hillas function that is to be fit
    x = (depth-Xo)/lambda_value
    m = (Xmax-Xo)/lambda_value
    n = ((x/m)**m) * np.exp(m-x)
    return np.log(n)

def get_chi2(x_true,x_predicted):
    #Returns the chi2 value
    x_predicted = np.array(x_predicted)
    x_true = np.array(x_true)
    chi2 = np.sum([((i-j)**2)/abs(j) for i,j in zip(x_true,x_predicted)])
    return chi2

def process_files(input_file):
    our_map = {
        'mass': [], 'energy': [], 'zenith': [], 'S125': [], 'energy_loss': [], 'he_stoch': [], 'he_stoch2': [],
        'Xmax': [], 'Xo': [], 'A': [], 'm_125': [], 'm_r': [], 'm_s': [], 'm_s2': [], 'm_o': [], 'm_chi2': [],
        'fit_status_m': [], 's_r': [], 's_o': [], 's_mean': [], 's_std': [], 's_chi2': [], 'fit_status_s': [],
        'charge': [], 'N': [], 'ghRedChiSqr': [], 'firstint': [], 'max_check': [], 'new_xmax': [], 'new_xo': [],
        'new_lambda': [], 'fit_status': [], 'new_chi2': [], 'difference': [], 'MaxNum': [], 'waveform_weight': []
    }

    l3_file = dataio.I3File(input_file, 'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics()

        if not np.allclose(l3_fr['IT73AnalysisInIceQualityCuts'].values(), 1):
            continue

        if 'Millipede_dEdX' in l3_fr:
            our_map['energy_loss'].append(l3_fr['Stoch_Reco'].eLoss_1500)
            our_map['he_stoch'].append(l3_fr['Stoch_Reco'].nHEstoch)
            our_map['he_stoch2'].append(l3_fr['Stoch_Reco2'].nHEstoch)
        else:
            our_map['energy_loss'].append(None)
            our_map['he_stoch'].append(None)
            our_map['he_stoch2'].append(None)

        our_map['mass'].append(4)
        our_map['energy'].append(l3_fr['MCPrimary'].energy)
        our_map['Xmax'].append(l3_fr['MCPrimaryInfo'].ghMaxDepth)
        our_map['Xo'].append(l3_fr['MCPrimaryInfo'].ghStartDepth)
        our_map['firstint'].append(l3_fr['MCPrimaryInfo'].firstIntDepth)
        our_map['ghRedChiSqr'].append(l3_fr['MCPrimaryInfo'].ghRedChiSqr)
        our_map['MaxNum'].append(l3_fr['MCPrimaryInfo'].ghMaxNum)

        try:
            our_map['waveform_weight'].append(l3_fr['IceTopWaveformWeight'].value)
        except KeyError:
            our_map['waveform_weight'].append(0)

        depth = []
        value = []

        for i in l3_fr['MCPrimaryInfo'].longProfile:
            depth.append(i.depth)
            value.append(i.numEMinus + i.numEPlus)

        our_map['max_check'].append(depth[np.argmax(value)] - l3_fr['MCPrimaryInfo'].ghMaxDepth)

        value = np.array(value)
        depth = np.array(depth)

        try:
            fit_params, _ = curve_fit(gh_function, xdata=depth[value > 0],
                                      ydata=np.log(value[value > 0] / max(value)), p0=[1, depth[value > 0][np.argmax(value[value > 0])], 2])

            our_map['new_xmax'].append(fit_params[1])
            our_map['new_xo'].append(fit_params[0])
            our_map['new_lambda'].append(fit_params[2])

            our_map['fit_status'].append(1 if np.isfinite(fit_params[0]) else 0)
            chi2_value = get_chi2(gh_function(depth[value > 0], fit_params[0], fit_params[1], fit_params[2]),
                                  np.log(value[value > 0] / max(value)))
            our_map['new_chi2'].append(chi2_value)
            our_map['difference'].append(fit_params[1] - depth[np.argmax(value)])

        except (ValueError, RuntimeError):
            our_map['new_xmax'].append(0)
            our_map['new_xo'].append(0)
            our_map['new_lambda'].append(0)
            our_map['fit_status'].append(0)
            our_map['new_chi2'].append(0)
            our_map['difference'].append(0)

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

        our_map['A'].append(l3_fr['CurvatureOnlyParams'].value(LaputopParameter.CurvParabA))

    l3_file.close()

    return our_map

count = 0

for i in my_list:
    print(i)
    if count == 0:
        new_map = process_files(i)
    else:
        event = process_files(i)
        for id,value in event.items():
            new_map[id]+=value 
    count+=1

df = pd.DataFrame(new_map)
df.to_csv(output_file)
