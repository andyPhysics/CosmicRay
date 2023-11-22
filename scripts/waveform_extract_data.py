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
from glob import glob
from shutil import copyfile

## Are you using tableio?                                                                                                                                                                 
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",type=str,default='IC86.2011',
                    dest="directory_number",help="directory number")
args = parser.parse_args()

#All files are in one directory for level3 files
l3_directory = '/data/ana/CosmicRay/IceTop_level3/exp/'
sub_dir1 = os.listdir(l3_directory+ '%s/'%(args.directory_number))

def process_files(my_input):
    print(my_input[0])
    outdir = '/data/user/amedina/CosmicRay/I3_files_data/%s'%(args.directory_number)+'/'+my_input[1]+'/' + my_input[2]

    assert os.path.isfile(my_input[0]), "file doesn't exist: {}".format(my_input[0])
    outfile = dataio.I3File(os.path.join(outdir,os.path.basename(my_input[0]).replace('.i3','.i3')),'w')
    l3_file = dataio.I3File(my_input[0],'r')

    while l3_file.more():
        l3_fr = l3_file.pop_physics()

        passed_all = True
        for cut in l3_fr['IT73AnalysisIceTopQualityCuts'].values():
            if not cut:
                passed_all = False

        if not passed_all:
            continue
        
        if 'IceTopWaveformWeight' in l3_fr:
            if l3_fr['IceTopWaveformWeight'].value==0:
                continue
        else:
            continue

        l3_header = l3_fr["I3EventHeader"]
        
        for l3_fr1 in l3_file.get_current_frame_and_deps():
            outfile.push(l3_fr1)
        

    outfile.close()
    l3_file.close()

my_list = []
os.mkdir('/data/user/amedina/CosmicRay/I3_files_data/%s/'%(args.directory_number))
for i in sub_dir1:
    if i == 'production':
        continue
    else:
        sub_dir2 = os.listdir(l3_directory+ '%s/'%(args.directory_number) + i)
        for j in sub_dir2:
            sub_dir3 = os.listdir(l3_directory+ '%s/'%(args.directory_number) + i + '/' + j)
            new_list = []
            for k in sub_dir3:
                if int(k[-1]) == 0:
                    new_list.append(k)
            if len(new_list)!= 0:
                try:
                    os.mkdir('/data/user/amedina/CosmicRay/I3_files_data/%s/'%(args.directory_number)+i)
                except FileExistsError:
                    pass

                try:
                    os.mkdir('/data/user/amedina/CosmicRay/I3_files_data/%s/'%(args.directory_number)+i+'/' + j)
                except FileExistsError:
                    pass
                for file_dir in new_list:
                    file_list = os.listdir(l3_directory+'%s/'%(args.directory_number) + i + '/' + j + '/' + file_dir)
                    for file in file_list:
                        if file.split('.')[-1] == 'gz':
                            file_name = l3_directory+'%s/'%(args.directory_number) + i + '/' + j + '/' + file_dir + '/' + file
                            copied_file = '/data/user/amedina/CosmicRay/I3_files_data/%s'%(args.directory_number)+'/'+i+'/' + j + '/' +file
                            copyfile(file_name,copied_file)
                        elif (file.split('.')[-1] == 'bz2') & (file.split('.')[-2] == 'i3'):
                            my_list.append([l3_directory+'%s/'%(args.directory_number) + i + '/' + j + '/' + file_dir + '/' + file,i,j])


pool = mp.Pool(processes=5)
pool.map(process_files,my_list)
