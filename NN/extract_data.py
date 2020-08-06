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

import numpy as np
import pandas as pd
from glob import glob

directory = '/data/user/amedina/CosmicRay/Analysis/'
number = sys.argv[1]
my_list = glob(directory + '/%s_*.i3.bz2'%(number))
output_file = sys.argv[2]
print(my_list)

def process_files(input_file):
    our_map = {}
    mass = []
    energy = []
    Xmax = []
    Xo = []

    zenith = []
    S125 = []
    energy_loss = []
    he_stoch = []
    he_stoch2 = []
    
    A = []

    m_125 = []
    m_z = []
    m_r = []
    m_o = []
    m_chi2 = []

    s_o = []
    s_1 = []
    s_mean = []
    s_std = []
    s_chi2 = []

    charge = []
    N = []
    
    l3_file = dataio.I3File(input_file,'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics() #This gets the frame
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

        zenith.append(l3_fr['Laputop'].dir.zenith)
        S125.append(l3_fr['LaputopParams'].value(LaputopParameter.Log10_S125))
        m_125.append(l3_fr['m_fit']['m_125'])
        m_o.append(l3_fr['m_fit']['m_o'])
        m_r.append(l3_fr['m_fit']['m_r'])
        m_z.append(l3_fr['m_fit']['m_z'])
        m_chi2.append(l3_fr['m_fit']['chi2'])

        s_o.append(l3_fr['s_fit']['s_o'])
        s_1.append(l3_fr['s_fit']['s_1'])
        s_mean.append(l3_fr['s_fit']['s_mean'])
        s_std.append(l3_fr['s_fit']['s_std'])
        s_chi2.append(l3_fr['s_fit']['chi2'])

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
        
        A.append(l3_fr['Laputop_newParams'].value(LaputopParameter.CurvParabA))
        
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
    our_map['m_z'] = m_z
    our_map['m_r'] = m_r
    our_map['m_o'] = m_o
    our_map['m_chi2'] = m_chi2
    our_map['s_o'] = s_o
    our_map['s_1'] = s_1
    our_map['s_mean'] = s_mean
    our_map['s_std'] = s_std
    our_map['s_chi2'] = s_chi2
    our_map['charge'] = charge
    our_map['N'] = N


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
