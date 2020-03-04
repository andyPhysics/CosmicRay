#!/home/amedina/build_stable/bin/python
import argparse
import os
import sys,getopt
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

from icecube.icetop_Level3_scripts.segments import level3_IceTop, level3_Coinc, IceTopQualityCuts
from icecube.icetop_Level3_scripts import icetop_globals



## Are you using tableio?
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

import numpy as np

## Set the log level
## Uncomment this in if you're doing code development!
icetray.set_log_level(icetray.I3LogLevel.LOG_INFO)
icetray.set_log_level_for_unit('Laputop', icetray.I3LogLevel.LOG_DEBUG)
icetray.set_log_level_for_unit('Curvature', icetray.I3LogLevel.LOG_DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",type=str,default='12360',
                    dest="directory_number",help="directory number")
args = parser.parse_args()


workspace = expandvars("$I3_BUILD")
#HOME = expandvars("$HOME")
HOME = '/data/user/amedina/Waveform_file'

data_set_number = args.directory_number
directory = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/'%(data_set_number)
directory_list = os.listdir('/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/'%(data_set_number))
file_list = []
for i in directory_list:
    x = np.array(os.listdir(directory+i))
    y = [directory+i+'/'+m for m in x]
    file_list.append(y)

file_list_all = np.hstack(file_list)
file_list_all = np.sort(file_list_all)[1:-1]
file_list_all = [str(i) for i in file_list_all]

file_list_all = np.array(np.array_split(file_list_all,1000))
    
#### PUT YOUR FAVORITE GCD AND INPUT FILE HERE

GCDfile = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/0000000-0000999/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz'%(data_set_number)
x_list = list(range(file_list_all.shape[0]))
x_list = np.array(x_list)

class Remove_WaveformRange(icetray.I3Module):
    def DAQ(self,frame):
        raw_data = True
        #if "CleanIceTopRawData" in frame:
        if not "CalibratedWaveformRange" in frame:
            raw_data = False
        if raw_data:
            del frame["CalibratedWaveformRange"]
            self.PushFrame(frame)
        else:
            self.PushFrame(frame)

    def Physics(self, frame):
        raw_data = True
        try:
            frame["CleanIceTopRawData"]
        except:
            raw_data = False
        if raw_data:
            del frame["CalibratedWaveformRange"]
            self.PushFrame(frame)
        else:
            self.PushFrame(frame)

@icetray.traysegment
def ExtractWaveforms(tray, name, InputLaunches=icetop_globals.icetop_raw_data,
                     OutputWaveforms='IceTopCalibratedWaveforms',
                     If=lambda f: True):
    OutputWaveforms = name + 'CalibratedWaveforms'
    OutputVEMWaveforms = name + 'VEMCalibratedWaveforms'

    tray.AddModule("I3WaveCalibrator", name+"_WaveCalibrator_IceTop",
                   Launches=InputLaunches, If=lambda fr: If(fr) and InputLaunches in fr,
                   Waveforms='ReextractedWaveforms',
                   WaveformRange="",
                   Errata="ReextractedErrata")

    tray.AddModule('I3WaveformSplitter', name + '_IceTopSplitter',
                   Input = 'ReextractedWaveforms',
                   HLC_ATWD = OutputWaveforms,
                   HLC_FADC = 'ReextractedHLCFADCWaveforms',
                   SLC = 'ReextractedSLCWaveforms',
                   Force = True, # ! put all maps in the frame, even if they are empty                                                                                                              
                   PickUnsaturatedATWD = True,  # ! do not keep all ATWDs, but only the highest non-saturated gain one                                                                              
                   If = lambda fr: If(fr) and InputLaunches in fr,
                   )

def output_i3_root(i):
    infile = np.array(file_list_all[i])


    first = infile[0].split('.')[-3]
    last = infile[-1].split('.')[-3]

    file_name = 'Subset_%s_%s'%(first,last)

    #### NAME YOUR OUTPUT FILES HERE

    I3_OUTFILE = HOME + '/' + data_set_number + '/' + file_name +'_' + data_set_number+ '.i3'
    ROOTFILE = HOME + '/' + data_set_number + '/' + file_name +'_' + data_set_number+ '.root'

    tray = I3Tray()

    ########## SERVICES FOR GULLIVER ##########

    datareadoutName="IceTopHLCSeedRTPulses"
    badtanksName= "IceTopHLCSeedRTExcludedTanks"

    #------------------- LET'S RUN SOME MODULES!  ------------------

    #**************************************************
    #                 Reader and whatnot
    #**************************************************
    infile=list(infile)

    tray.AddModule("I3Reader","reader")(
        ("FileNameList", [GCDfile]+infile)
    )

    tray.AddModule(Remove_WaveformRange)

    from icecube.icetop_Level3_scripts.functions import count_stations
    from icecube.icetop_Level3_scripts.modules import FilterWaveforms

    tray.AddModule(FilterWaveforms, 'FilterWaveforms',   #Puts IceTopWaveformWeight in the frame.                                                     
                   pulses=icetop_globals.icetop_hlc_pulses,
                   If = lambda frame: icetop_globals.icetop_hlc_pulses in frame and count_stations(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, icetop_globals.icetop_hlc_pulses)) >= 5)

    tray.AddSegment(ExtractWaveforms, 'IceTop')
                       #If= lambda frame: "IceTopWaveformWeight" in frame and frame["IceTopWaveformWeight"].value!=0)

    # Extract HLC pulses
    tray.AddModule('I3TopHLCPulseExtractor', 'TopHLCPulseExtractor',
                   PEPulses  = 'IceTopHLCPEPulses',         # Pulses in PE, set to empty string to disable output
                   PulseInfo = 'IceTopHLCPulseInfo',        # PulseInfo: amplitude, rise time, etc. Empty string to disable
                   VEMPulses = 'IceTopHLCVEMPulses',        # Pulses in VEM, set to empty string to disable
                   Waveforms = 'IceTopCalibratedWaveforms',   # Input HLC waveforms from WaveCalibrator
                   BadDomList = "BadDomsList"
               )

    # Extract SLC pulses
    tray.AddModule('I3TopSLCPulseExtractor', 'TopSLCPulseExtractor',
                   PEPulses  = 'IceTopSLCPEPulses',         # (see above ...)
                   VEMPulses = 'IceTopSLCVEMPulses',
                   Waveforms = 'ReextractedSLCWaveforms',   # Input SLC waveforms from WaveCalibrator
                   BadDomList = "BadDomsListSLC"
               )


    itpulses='IceTopHLCSeedRTPulses'

    #Keys to Keep

    wanted_general=['I3EventHeader',
                    icetop_globals.filtermask,
                    'I3TriggerHierarchy']

    wanted_MC = ['I3MCPESeriesMap',
                 'I3MCTree'] 


    wanted_icetop_waveforms=['IceTopCalibratedWaveforms',
                             'IceTopWaveformWeight',
                             'ReextractedHLCFADCWaveforms',
                             'ReextractedSLCWaveforms',
                             'IceTopHLCPEPulses',
                             'IceTopHLCPulseInfo',
                             'IceTopHLCVEMPulses',
                             'IceTopSLCPEPulses',                                                                                                                      
                             'IceTopSLCVEMPulses',
                             'IceTopComponentPulses_Electron',
                             'IceTopComponentPulses_ElectronFromChargedMesons', 
                             'IceTopComponentPulses_Gamma',
                             'IceTopComponentPulses_GammaFromChargedMesons',
                             'IceTopComponentPulses_Hadron',  
                             'IceTopComponentPulses_Muon']

    wanted=wanted_general+wanted_icetop_waveforms
    
    tray.AddModule("Keep", 'DropObjects',
                   Keys = wanted
                   )

    
    tray.AddModule("I3Writer","EventWriter")(
        ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
        ("Filename",I3_OUTFILE),
    )

    ## Output root file
    root = I3ROOTTableService(ROOTFILE,"aTree")
    tray.AddModule(I3TableWriter,'writer')(
        ("tableservice", root),
        ("keys",['I3EventHeader',
                 'IceTopCalibratedWaveforms',
                 'IceTopWaveformWeight',
                 'ReextractedHLCFADCWaveforms',
                 'ReextractedSLCWaveforms',
                 'IceTopHLCPEPulses',
                 'IceTopHLCPulseInfo',
                 'IceTopHLCVEMPulses',
                 'IceTopSLCPEPulses',
                 'IceTopSLCVEMPulses',
                 'I3MCPESeriesMap',
                 'I3MCTree',
                 'IceTopComponentPulses_Electron',
                 'IceTopComponentPulses_ElectronFromChargedMesons',
                 'IceTopComponentPulses_Gamma',
                 'IceTopComponentPulses_GammaFromChargedMesons',
                 'IceTopComponentPulses_Hadron',
                 'IceTopComponentPulses_Muon'
             ]),
        ("subeventstreams",['InIceSplit',"ice_top"])
    )

   
    # Execute the Tray
    # Just to make sure it's working!
    tray.Execute()

pool = mp.Pool(processes=5)
pool.map(output_i3_root,x_list)
