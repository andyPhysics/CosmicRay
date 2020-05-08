#!/home/amedina/build_stable/bin/python                                                                                                                                                   
import argparse
import os,sys,getopt

from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses,simclasses,tpx

from icecube.icetray import I3Module, OMKey
from icecube.dataclasses import I3EventHeader, I3Particle
from icecube.recclasses import I3LaputopParams

import numpy as np
import pandas as pd


data_set_number = str(12360)
filenumber = str(0)
file_name = '/data/user/amedina/CosmicRay/Analysis/%s_%s.i3.bz2'%(data_set_number,filenumber)
l3_file = dataio.I3File(file_name,'r')

event = {}

count = 0

while l3_file.more():
    frame = l3_file.pop_physics()
    if (count > 100):
        break
    elif (np.log10(frame['MCPrimary'].energy) < 7):
        continue

    eventinfo = {}
    for omkey in frame['WaveformInfo'].keys():
        eventinfo[omkey] = {}
        eventinfo[omkey]['run_number'] = frame['I3EventHeader'].run_id
        eventinfo[omkey]['event_number'] = frame['I3EventHeader'].event_id
        eventinfo[omkey]['zenith'] = frame['MCPrimary'].dir.zenith
        eventinfo[omkey]['azimuth'] = frame['MCPrimary'].dir.azimuth
        eventinfo[omkey]['energy'] = np.log10(frame['MCPrimary'].energy)
        dom = int(omkey.split(',')[1])
        string = int(omkey.split(',')[0].split('(')[1])
        position = frame['I3Geometry'].omgeo[OMKey(string,dom)].position
        eventinfo[omkey]['x'] = position.x
        eventinfo[omkey]['y'] = position.y
        eventinfo[omkey]['z'] = position.z
        eventinfo[omkey]['m'] = frame['WaveformInfo'][omkey]['m']
        eventinfo[omkey]['s'] = frame['WaveformInfo'][omkey]['s']
        eventinfo[omkey]['t0'] = frame['WaveformInfo'][omkey]['t_0']
        eventinfo[omkey]['charge'] = frame['WaveformInfo'][omkey]['Charge']
        eventinfo[omkey]['chargePe'] = frame['WaveformInfo'][omkey]['Charge_PE']
        eventinfo[omkey]['chargeVEM'] = frame['LaputopHLCVEM'][OMKey(string,dom)][0].charge
        eventinfo[omkey]['chi2'] = frame['WaveformInfo'][omkey]['chi2']
        eventinfo[omkey]['sigmat0'] = frame['WaveformInfo'][omkey]['sigma_t']
        eventinfo[omkey]['sigmam'] = frame['WaveformInfo'][omkey]['sigma_s']
        eventinfo[omkey]['sigmas'] = frame['WaveformInfo'][omkey]['sigma_m']
        fe_impedance = frame['I3Calibration'].dom_cal[OMKey(string,dom)].front_end_impedance
        spe_mean = dataclasses.spe_mean(frame['I3DetectorStatus'].dom_status[OMKey(string,dom)],frame['I3Calibration'].dom_cal[OMKey(string,dom)])
        eventinfo[omkey]['feimpedance'] = fe_impedance
        eventinfo[omkey]['spemean'] = spe_mean
        pe_per_vem = frame['I3Calibration'].vem_cal[OMKey(string,dom)]
        eventinfo[omkey]['pe_per_vem'] = pe_per_vem.pe_per_vem/pe_per_vem.corr_factor


    event['event_%s'%(count)] = eventinfo
        
    count+=1

l3_file.close()

user_ids = []
frames = []

for user_id, d in event.items():
    user_ids.append(user_id)
    if 'event' in user_id:
        frames.append(pd.DataFrame.from_dict(d,orient='index'))
    

df = pd.concat(frames, keys=user_ids)

df.to_csv('Events.csv')
print(df.head())
