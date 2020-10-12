#!/home/amedina/build_stable/bin/python                                                                                                                                                   
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
    mass = []
    energy = []
    Xmax = []
    MaxNum = []
    Xo = []

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
    ghredchi = []
    firstint = []
    max_check = []

    new_xmax = []
    new_xo = []
    new_lambda = []
    fit_status = []
    chi2 = []
    difference = []

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

        

        mass.append(4)
        energy.append(l3_fr['MCPrimary'].energy)
        Xmax.append(l3_fr['MCPrimaryInfo'].ghMaxDepth)
        Xo.append(l3_fr['MCPrimaryInfo'].ghStartDepth)
        firstint.append(l3_fr['MCPrimaryInfo'].firstIntDepth)
        ghredchi.append(l3_fr['MCPrimaryInfo'].ghRedChiSqr)
        MaxNum.append(l3_fr['MCPrimaryInfo'].ghMaxNum)
        waveform_weight.append(l3_fr['IceTopWaveformWeight'].value)
        depth = []
        value = []
        for i in l3_fr['MCPrimaryInfo'].longProfile:
            depth.append(i.depth)
            value.append(i.numEMinus+i.numEPlus)

        max_check.append(depth[np.argmax(value)]-l3_fr['MCPrimaryInfo'].ghMaxDepth)


        value = np.array(value)
        depth = np.array(depth)
        try:
            fit = curve_fit(gh_function,xdata=depth[value>0],ydata=np.log(value[value>0]/max(value)),p0 = [1,depth[value>0][np.argmax(value[value>0])],2])
            new_xmax.append(fit[0][1])
            new_xo.append(fit[0][0])
            new_lambda.append(fit[0][2])
            if np.isfinite(fit[1][0][0]):
                fit_status.append(1)
            else:
                fit_status.append(0)
            chi2_value = get_chi2(gh_function(depth[value>0],fit[0][0],fit[0][1],fit[0][2]),np.log(value[value>0]/max(value)))
            chi2.append(chi2_value)
            difference.append(fit[0][1]-depth[np.argmax(value)])

        except (ValueError,RuntimeError) as err:
            new_xmax.append(0)
            new_xo.append(0)
            new_lambda.append(0)
            fit_status.append(0)
            chi2.append(0)
            difference.append(0)


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
    our_map['mass'] = mass
    our_map['energy'] = energy
    our_map['zenith'] = zenith
    our_map['S125'] = S125
    our_map['energy_loss'] = energy_loss
    our_map['he_stoch'] = he_stoch
    our_map['he_stoch2'] = he_stoch2
    our_map['Xmax'] = Xmax
    our_map['Xo'] = Xo
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
    our_map['fit_status_s'] = fit_status_s
    our_map['charge'] = charge
    our_map['N'] = N
    our_map['ghRedChiSqr'] = ghredchi
    our_map['firstint'] = firstint
    our_map['max_check'] = max_check
    our_map['new_xmax'] = new_xmax
    our_map['new_xo'] = new_xo
    our_map['new_lambda'] = new_lambda
    our_map['fit_status'] = fit_status
    our_map['new_chi2'] = chi2
    our_map['difference'] = difference
    our_map['MaxNum'] = MaxNum
    our_map['waveform_weight'] = waveform_weight


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
