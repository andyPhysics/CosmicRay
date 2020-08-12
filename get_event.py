#!/home/amedina/build_stable/bin/python                                                                                                                                                   
import argparse
import os,sys,getopt

from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses,simclasses,tpx

from icecube.icetray import I3Module, OMKey
from icecube.dataclasses import I3EventHeader, I3Particle
from icecube.recclasses import I3LaputopParams, LaputopParameter

import numpy as np
import pandas as pd


data_set_number = str(12360)
filenumber = str(0)
file_name = '/data/user/amedina/CosmicRay/Analysis/%s_%s.i3.bz2'%(data_set_number,filenumber)
l3_file = dataio.I3File(file_name,'r')
output_name = 'Events.csv'

event = {}

count = 0

while l3_file.more():
    frame = l3_file.pop_physics()
    charge = []
    time = []
    for i in frame['LaputopHLCVEM'].keys():
        charge.append(frame['LaputopHLCVEM'][i][0].charge)
    try:
        check = np.array(charge)>10
    except RuntimeWarning:
        continue
    if (count > 1000):
        break
    elif (np.log10(frame['MCPrimary'].energy) < 7)&(np.sum(check)>5):
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
        eventinfo[omkey]['ShowerCOG_x'] = frame['Laputop'].pos.x
        eventinfo[omkey]['ShowerCOG_y'] = frame['Laputop'].pos.y
        eventinfo[omkey]['ShowerCOG_z'] = frame['Laputop'].pos.z
        eventinfo[omkey]['ShowerCOG_time'] = frame['ShowerCOG'].time
        eventinfo[omkey]['ShowerCOG_zen'] = frame['ShowerPlane'].dir.zenith
        eventinfo[omkey]['ShowerCOG_az'] = frame['ShowerPlane'].dir.azimuth
        eventinfo[omkey]['m'] = frame['WaveformInfo'][omkey]['m']
        eventinfo[omkey]['s'] = frame['WaveformInfo'][omkey]['s']
        eventinfo[omkey]['charge'] = frame['WaveformInfo'][omkey]['Charge']
        eventinfo[omkey]['chargePe'] = frame['WaveformInfo'][omkey]['Charge_PE']
        eventinfo[omkey]['chargeVEM'] = frame['LaputopHLCVEM'][OMKey(string,dom)][0].charge
        eventinfo[omkey]['chi2'] = frame['WaveformInfo'][omkey]['chi2']
        eventinfo[omkey]['sigmam'] = frame['WaveformInfo'][omkey]['sigma_s']
        eventinfo[omkey]['sigmas'] = frame['WaveformInfo'][omkey]['sigma_m']
        eventinfo[omkey]['A'] = frame['Laputop_newParams'].value(LaputopParameter.CurvParabA)
        fe_impedance = frame['I3Calibration'].dom_cal[OMKey(string,dom)].front_end_impedance
        spe_mean = dataclasses.spe_mean(frame['I3DetectorStatus'].dom_status[OMKey(string,dom)],frame['I3Calibration'].dom_cal[OMKey(string,dom)])
        eventinfo[omkey]['feimpedance'] = fe_impedance
        eventinfo[omkey]['spemean'] = spe_mean
        pe_per_vem = frame['I3Calibration'].vem_cal[OMKey(string,dom)]
        eventinfo[omkey]['pe_per_vem'] = pe_per_vem.pe_per_vem/pe_per_vem.corr_factor
        eventinfo[omkey]['angular_resolution'] = frame['LaputopParams'].angular_resolution
        eventinfo[omkey]['chi2_ldf'] = frame['LaputopParams'].chi2_ldf
        eventinfo[omkey]['chi2_time'] = frame['LaputopParams'].chi2_time
        eventinfo[omkey]['Laputop_dir_zenith'] = frame['Laputop'].dir.zenith
        eventinfo[omkey]['Laputop_dir_azimuth'] = frame['Laputop'].dir.azimuth
        eventinfo[omkey]['Laputop_new_zenith'] = frame['Laputop_new'].dir.zenith
        eventinfo[omkey]['Laputop_new_azimuth'] = frame['Laputop_new'].dir.azimuth
        eventinfo[omkey]['Laputop_time'] = frame['Laputop'].time
        eventinfo[omkey]['Laputop_pos_x'] = frame['Laputop'].pos.x
        eventinfo[omkey]['Laputop_pos_y'] = frame['Laputop'].pos.y
        eventinfo[omkey]['t_0'] = frame['LaputopHLCVEM'][OMKey(string,dom)][0].time
        eventinfo[omkey]['Xmax'] = frame['MCPrimaryInfo'].ghMaxDepth



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

df.to_csv(output_name)
print(df.head())
