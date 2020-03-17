#!/home/amedina/build_stable/bin/python
import argparse
import os
import sys,getopt
#import logging
from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, frame_object_diff,simclasses,WaveCalibrator,tpx
load('millipede')
load('stochastics')

from icecube.icetray import I3Module
from icecube.dataclasses import I3EventHeader, I3Particle
from icecube.recclasses import I3LaputopParams
from icecube.stochastics import *
from icecube.tpx import *
import multiprocessing as mp

## Are you using tableio?
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService

import numpy as np

## Set the log level
## Uncomment this in if you're doing code development!
icetray.set_log_level(icetray.I3LogLevel.LOG_INFO)
icetray.set_log_level_for_unit('Laputop', icetray.I3LogLevel.LOG_DEBUG)
icetray.set_log_level_for_unit('Curvature', icetray.I3LogLevel.LOG_DEBUG)

c = .299 #m/ns

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

I3_OUTFILE = output_directory + data_set_number + '.i3.bz2' 
ROOTFILE = output_directory + data_set_number + '.root'
    
#### PUT YOUR FAVORITE GCD AND INPUT FILE HERE

GCDFile = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/0000000-0000999/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz'%(data_set_number)


@icetray.traysegment
def ExtractWaveforms(tray, name, InputLaunches='CleanIceTopRawData',
                     OutputWaveforms='CalibratedHLCWaveforms'):

    tray.AddModule("I3WaveCalibrator", name+"_WaveCalibrator_IceTop",
                   Launches=InputLaunches,
                   Waveforms='ReextractedWaveforms',
                   Errata="ReextractedErrata")

    tray.AddModule('I3WaveformSplitter', name + '_IceTopSplitter',
                   Input = 'ReextractedWaveforms',
                   HLC_ATWD = OutputWaveforms,
                   SLC = 'CalibratedSLCWaveforms',
                   PickUnsaturatedATWD = True,  # ! do not keep all ATWDs, but only the highest non-saturated gain one                                        
                   )

def GetRiseTime(values,time,percentage):
    bin_value = 0
    for i in zip(values,time):
        if i[0] > max(values) * percentage:
            return i[1],bin_value
        bin_value+=1

class Process_Waveforms(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Physics(self, frame):
        mask = frame['IceTopLaputopSeededSelectedHLC'] #These are the keys of the HLC events that are used for reconstruction
        check = False
        if 'IceTopLaputopSeededSelectedSLC' in frame:
            check = True
            mask_2 = frame['IceTopLaputopSeededSelectedSLC']
        waveforms = frame['CalibratedHLCWaveforms']
        keys = []
        keys_2 = []
        for i in zip(np.hstack(mask.bits),frame[mask.source].keys()):
            if i[0]:
                keys.append(i[1])
        
        if check and len(mask_2.bits) != 0:
            for i in zip(np.hstack(mask_2.bits),frame[mask_2.source].keys()):
                if i[0]:
                    keys_2.append(i[1])

        HLCWaveforms = dataclasses.I3WaveformSeriesMap()
        PE = dataclasses.I3RecoPulseSeriesMap()
        PE_2 = dataclasses.I3RecoPulseSeriesMap()

        for i in keys:
            HLCWaveforms[i] = waveforms[i]
            PE[i] = frame['IceTopHLCPEPulses'][i]

        if check:
            for i in keys_2:
                PE_2[i] = frame['IceTopSLCPEPulses'][i]
            
        frame['LaputopHLCWaveforms'] = HLCWaveforms
        frame['LaputopHLCPE'] = PE
        frame['LaputopSLCPE'] = PE_2
        self.PushFrame(frame)

class Extract_info(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Physics(self, frame):
        OutputWaveformInfo = dataclasses.I3MapStringStringDouble()
        
        for i in frame['LaputopHLCWaveforms'].keys():
            key = str(i)
            waveform = frame['LaputopHLCWaveforms'][i][0].waveform
            time = frame['LaputopHLCWaveforms'][i][0].time
            binwidth = frame['LaputopHLCWaveforms'][i][0].binwidth
            
            time_values = [i*binwidth + time for i in range(len(waveform))]
            waveform_check = np.array(waveform) >= 0

            CDF = []
            CDF_value = 0
            for j in range(len(waveform)):
                if waveform_check[j]:
                    CDF_value = CDF_value + binwidth * waveform[j]
                    CDF.append(CDF_value)
                else:
                    CDF.append(CDF_value)


            fe_impedance = frame['I3Calibration'].dom_cal[i].front_end_impedance
            charge = CDF[-1]/fe_impedance
            spe_mean = dataclasses.spe_mean(frame['I3DetectorStatus'].dom_status[i],frame['I3Calibration'].dom_cal[i])

            charge_pe = charge/spe_mean


            Amplitude = 0
            for k in waveform:
                if k-Amplitude >= 0:
                    Amplitude = k
                else:
                    break

            Time_10,bin_10 = GetRiseTime(CDF,time_values,0.1)
            Time_50,bin_50 = GetRiseTime(CDF,time_values,0.5)
            Time_90,bin_90 = GetRiseTime(CDF,time_values,0.9)

            ninety_slope = (CDF[bin_90] - CDF[bin_10])/(bin_90 - bin_10)
            fifty_slope = (CDF[bin_50] - CDF[bin_10])/(bin_50 - bin_10)

            tmin = -200.0 #ns

            leading_edge = time + (bin_10 - CDF[bin_10]/ninety_slope)*binwidth
            if ninety_slope <= 0 or not np.isfinite(ninety_slope) or not (leading_edge >= time + tmin):
                leading_edge = time + tmin
                                                            

            OutputWaveformInfo[key] = dataclasses.I3MapStringDouble()
            OutputWaveformInfo[key]['StartTime'] = time
            OutputWaveformInfo[key]['Binwidth'] = binwidth
            OutputWaveformInfo[key]['Time_50'] = Time_50
            OutputWaveformInfo[key]['Time_90'] = Time_90
            OutputWaveformInfo[key]['Time_10'] = Time_10
            OutputWaveformInfo[key]['Amplitude'] = Amplitude
            OutputWaveformInfo[key]['Charge'] = charge
            OutputWaveformInfo[key]['Charge_PE'] = charge_pe
            OutputWaveformInfo[key]['90_slope'] = ninety_slope
            OutputWaveformInfo[key]['leading_edge'] = leading_edge

        frame['WaveformInfo'] = OutputWaveformInfo
        self.PushFrame(frame)


class Get_data(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Physics(self, frame):
        a = 4.823 * 10.0**(-4.0) #ns/m^2
        b = 19.41 #ns
        sigma = 83.5 #m

        output_map = dataclasses.I3RecoPulseSeriesMap()

        Laputop = frame['Laputop']
        x_core = np.array([Laputop.pos.x,Laputop.pos.y,Laputop.pos.z])
        t_core = Laputop.time
        theta = Laputop.dir.theta 
        phi = Laputop.dir.phi
        n = np.array([np.sin(theta) * np.cos(phi) , np.sin(theta) * np.sin(phi), np.cos(theta)])
        
        for i in frame['LaputopSLCPE'].keys():
            output_map[i] = dataclasses.I3RecoPulseSeries()

            pulse = dataclasses.I3RecoPulse()
            pulse.charge = frame['LaputopSLCPE'][i][0].charge
            time = frame['LaputopSLCPE'][i][0].time
            position_dom = frame['I3Geometry'].omgeo[i].position
            x_dom = np.array([position_dom.x , position_dom.y , position_dom.z])
            
            R_square = np.dot(x_core-x_dom,x_core-x_dom)
            delta_T = a * R_square + b * (1 - np.exp(-R_square/(2*(sigma**2.0))))
            time_signal = time + (1/c) * np.dot(x_core-x_dom,n) + delta_T
            
            pulse.time = time_signal
            output_map[i].append(pulse)

        for i in frame['LaputopHLCPE'].keys():
            output_map[i] = dataclasses.I3RecoPulseSeries()
            pulse = dataclasses.I3RecoPulse()

            pulse.charge = frame['LaputopHLCPE'][i][0].charge
            time = frame['LaputopHLCPE'][i][0].time
            position_dom = frame['I3Geometry'].omgeo[i].position
            x_dom = np.array([position_dom.x , position_dom.y , position_dom.z])

            R_square = np.dot(x_core-x_dom,x_core-x_dom)
            delta_T = a * R_square + b * (1 - np.exp(-R_square/(2*(sigma**2.0))))
            time_signal = time + (1/c) * np.dot(x_core-x_dom,n) + delta_T

            pulse.time = time_signal
            output_map[i].append(pulse)

        frame['All_pulses'] = output_map
        self.PushFrame(frame)

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
               Waveforms = 'CalibratedHLCWaveforms',   # Input HLC waveforms from WaveCalibrator
               BadDomList = "BadDomsList"
)

# Extract SLC pulses
tray.AddModule('I3TopSLCPulseExtractor', 'TopSLCPulseExtractor',
               PEPulses  = 'IceTopSLCPEPulses',         # (see above ...)
               VEMPulses = 'IceTopSLCVEMPulses',
               Waveforms = 'CalibratedSLCWaveforms',   # Input SLC waveforms from WaveCalibrator
               BadDomList = "BadDomsListSLC"
)

tray.AddModule(Process_Waveforms,'Process_wavefomrs')

tray.AddModule(Extract_info)

tray.AddModule(Get_data)

tray.AddModule("I3Writer","EventWriter")(
    ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
    ("Filename",I3_OUTFILE),
)

## Output root file
root = I3ROOTTableService(ROOTFILE,"aTree")
tray.AddModule(I3TableWriter,'writer')(
    ("tableservice", root),
    ("keys",['I3EventHeader',
             'CalibratedHLCWaveforms',
             'CalibratedSLCWaveforms',
             'MCPrimary',
             'MCPrimaryInfo',
             'LaputopHLCWaveforms',
             'IceTopHLCPEPulses',
             'IceTopHLCPulseInfo',
             'IceTopHLCVEMPulses',
             'IceTopSLCPEPulses',
             'IceTopSLCVEMPulses',
             'LaputopHLCPE',
             'LaputopSLCPE',
             'LaputopParams',
             'All_data'
         ]),
    ("subeventstreams",['InIceSplit',"ice_top"])
)

   
# Execute the Tray
# Just to make sure it's working!
tray.Execute()

