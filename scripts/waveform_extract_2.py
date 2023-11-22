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
from icecube.recclasses import I3LaputopParams
from icecube.stochastics import *
import multiprocessing as mp

## Are you using tableio?                                                                                                                                                                 
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",type=str,default='12360',
                    dest="directory_number",help="directory number")
args = parser.parse_args()

#level2 files are separated into 1000 per sub directory
l2_directory = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/'%(args.directory_number) 

#All files are in one directory for level3 files
l3_directory = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/%s/'%(args.directory_number)

l3_file_list = []
l3_directory_list = os.listdir(l3_directory)
for i in l3_directory_list:
    if i.endswith('.i3.gz'):
        l3_file_list.append(l3_directory+i)

l3_file_list = np.array(l3_file_list)
l3_file_list = np.sort(l3_file_list)

check_numbers = [i.split('.')[-3].split('Run')[-1] for i in l3_file_list]

l2_directory_list = os.listdir(l2_directory)
l2_file_list = []
for i in l2_directory_list:
    x = np.array(os.listdir(l2_directory+i))
    y = [l2_directory+i+'/'+m for m in x]
    for j in y:
        if j.endswith('.i3.bz2'):
            if j.split('.')[-3] in check_numbers:
                l2_file_list.append(j)

l2_file_list = np.sort(l2_file_list)

def process_files(file_tuple):
    print(file_tuple)
    outdir = '/data/user/amedina/CosmicRay/I3_Files/%s'%(args.directory_number)

    assert os.path.isfile(file_tuple[1]), "file doesn't exist: {}".format(file_tuple[1])
    assert os.path.isfile(file_tuple[0]), "file doesn't exist: {}".format(file_tuple[0])
    outfile = dataio.I3File(os.path.join(outdir,os.path.basename(file_tuple[1]).replace('.i3','.i3')),'w')
    l3_file = dataio.I3File(file_tuple[1],'r')
    l2_file = dataio.I3File(file_tuple[0],'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics()

        passed_all = True
        for cut in l3_fr['IT73AnalysisIceTopQualityCuts'].values():
            if not cut:
                passed_all = False

        if not passed_all:
            continue

        l3_header = l3_fr["I3EventHeader"]
        while l2_file.more():
            try:
                l2_fr = l2_file.pop_daq()
            except RuntimeError:
                continue
            l2_header = l2_fr["I3EventHeader"]
            if l2_header.run_id == l3_header.run_id and l2_header.event_id == l3_header.event_id:
                for l3_fr1 in l3_file.get_current_frame_and_deps():
                    l3_fr1['CleanIceTopRawData'] = l2_fr['CleanIceTopRawData']
                    outfile.push(l3_fr1)
                break

    outfile.close()
    l3_file.close()
    l2_file.close()

file_list = list(zip(l2_file_list,l3_file_list))

pool = mp.Pool(processes=5)
pool.map(process_files,file_list)
