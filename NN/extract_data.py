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
my_list = glob(directory + '/12362_*.i3.bz2')
output_file = 'Iron.csv'
print(my_list)

def process_files(input_file):
    our_map = {}
    mass = []
    energy = []

    zenith = []
    S125 = []
    energy_loss = []
    he_stoch = []
    he_stoch2 = []
    
    alpha = []
    beta = []
    chi = []
    omega = []
    
    
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

        zenith.append(l3_fr['Laputop'].dir.zenith)
        S125.append(l3_fr['LaputopParams'].value(LaputopParameter.Log10_S125))
        
        #my_variables, beta is not related to the age of the shower
        alpha.append(l3_fr['my_fit']['alpha'])
        beta.append(l3_fr['my_fit']['beta'])
        chi.append(l3_fr['my_fit']['chi2'])
        omega.append(l3_fr['my_fit']['omega'])
        
    l3_file.close()
    our_map['mass'] = mass
    our_map['energy'] = energy
    our_map['zenith'] = zenith
    our_map['S125'] = S125
    our_map['energy_loss'] = energy_loss
    our_map['he_stoch'] = he_stoch
    our_map['he_stoch2'] = he_stoch2
    our_map['alpha'] = alpha
    our_map['beta'] = beta
    our_map['chi'] = chi
    our_map['omega'] = omega

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
