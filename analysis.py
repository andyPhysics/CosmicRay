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

output_directory = '/data/user/amedina/CosmicRay/Analysis/'

data_set_number = args.directory_number

directory = '/data/user/amedina/CosmicRay/I3_Files/%s/'%(data_set_number)
directory_list = os.listdir(directory)
file_list = []
for i in directory_list:
    file_list.append(directory+i)

file_list = np.array(file_list)

I3_OUTFILE = output_directory + data_set_number + '.i3.bz2' 
ROOTFILE = output_directory + data_set_number + '.root'
    
#### PUT YOUR FAVORITE GCD AND INPUT FILE HERE

GCDFile = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/0000000-0000999/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz'%(data_set_number)

@icetray.traysegment
def ExtractWaveforms(tray, name, InputLaunches='CleanIceTopRawData',
                     OutputWaveforms='IceTopCalibratedWaveforms'):
    OutputWaveforms = name + 'CalibratedWaveforms'

    tray.AddModule("I3WaveCalibrator", name+"_WaveCalibrator_IceTop",
                   Launches=InputLaunches,
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
                   )

class Process_Waveforms(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Physics(self, frame):
        mask = frame['IceTopLaputopSeededSelectedHLC'] #These are the keys of the HLC events that are used for reconstruction
        waveforms = frame['IceTopCalibratedWaveforms']
        keys = []
        for i in zip(np.hstack(mask.bits),frame[mask.source].keys()):
            if i[0]:
                keys.append(i[1])
        HLCWaveforms = dataclasses.I3WaveformSeriesMap()
        for i in keys:
            HLCWaveforms[i] = waveforms[i]
        frame['HLCWaveforms'] = HLCWaveforms
        self.PushFrame(frame)

file_list = ['/data/user/amedina/CosmicRay/I3_Files/12360/Level3_IC86.2012_12360_Run003010.i3.gz']

tray = I3Tray()

########## SERVICES FOR GULLIVER ##########

#------------------- LET'S RUN SOME MODULES!  ------------------

#**************************************************
#                 Reader and whatnot
#**************************************************

tray.AddModule("I3Reader","reader")(("FileNameList", [GCDFile] + file_list))

tray.AddSegment(ExtractWaveforms, 'IceTop')
               

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

tray.AddModule(Process_Waveforms)

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

