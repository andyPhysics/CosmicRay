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
        PE = frame['IceTopLaputopSeededSelectedHLC'].apply(frame)
        PE_2 = frame['IceTopLaputopSeededSelectedSLC'].apply(frame)

        waveforms = frame['CalibratedHLCWaveforms']

        HLCWaveforms = dataclasses.I3WaveformSeriesMap()

        for i in PE.keys():
            HLCWaveforms[i] = waveforms[i]
            

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
        t_cog = frame['ShowerCOG'].time
        zenith_core = frame['Laputop'].dir.zenith
        azimuth_core = frame['Laputop'].dir.azimuth
        unit_core = np.array([np.sin(zenith_core)*np.cos(azimuth_core),
                              np.sin(zenith_core)*np.sin(azimuth_core),
                              np.cos(zenith_core)])

        output_map = dataclasses.I3RecoPulseSeriesMap()
        output_10 = dataclasses.I3RecoPulseSeriesMap()

        Laputop = frame['Laputop']
        x_core = np.array([Laputop.pos.x,Laputop.pos.y,Laputop.pos.z])
        t_core = Laputop.time
        theta = Laputop.dir.theta 
        phi = Laputop.dir.phi
        n = np.array([np.sin(theta) * np.cos(phi) , np.sin(theta) * np.sin(phi), np.cos(theta)])

        radius = dataclasses.I3MapKeyVectorDouble()
        radius_old = dataclasses.I3MapKeyVectorDouble()

        for i in frame['LaputopSLCPE'].keys():
            output_map[i] = dataclasses.I3RecoPulseSeries()
            output_10[i] = dataclasses.I3RecoPulseSeries()

            pulse = dataclasses.I3RecoPulse()
            pulse2 = dataclasses.I3RecoPulse()

            vec = []
            vec_old = []
            
            for j in frame['LaputopSLCPE'][i]:
                pulse.charge = j.charge
                pulse2.charge = j.charge

                time = j.time
                position_dom = frame['I3Geometry'].omgeo[i].position
                x_dom = np.array([position_dom.x , position_dom.y , position_dom.z])
                
                Radius = np.dot(x_dom-x_core,x_dom-x_core)**0.5
                unit_dom = (x_dom-x_core)/Radius
                true_radius = np.dot(unit_dom-unit_core,x_dom-x_core)
                #time_signal = time - (1/c) * np.dot(x_core-x_dom,n) - t_cog
                time_signal = time
                vec.append(true_radius)
                vec_old.append(Radius)
                pulse.time = time_signal
                pulse2.time = time_signal
                output_map[i].append(pulse)
                output_10[i].append(pulse2)

            radius[i] = np.array(vec)
            radius_old[i] = np.array(vec_old)

        for i in frame['LaputopHLCPE'].keys():
            output_map[i] = dataclasses.I3RecoPulseSeries()
            output_10[i] = dataclasses.I3RecoPulseSeries()

            pulse = dataclasses.I3RecoPulse()
            pulse2 = dataclasses.I3RecoPulse()
            
            vec = []
            vec_old = []
            for j in frame['LaputopHLCPE'][i]:
                pulse.charge = j.charge
                pulse2.charge = j.charge
                time = j.time
                position_dom = frame['I3Geometry'].omgeo[i].position
                x_dom = np.array([position_dom.x , position_dom.y , position_dom.z])

                Radius = np.dot(x_dom-x_core,x_dom-x_core)**0.5
                unit_dom = (x_dom-x_core)/Radius
                true_radius = np.dot(unit_dom-unit_core,x_dom-x_core)

                #time_signal = time - (1/c) * np.dot(x_core-x_dom,n) - t_cog
                time_signal = time
            
                vec.append(true_radius)
                vec_old.append(Radius)
                pulse.time = time_signal
                key = str(i)
                pulse2.time = frame['WaveformInfo'][key]['Time_50'] - frame['WaveformInfo'][key]['Time_10'] + (1/c) * np.dot(x_core-x_dom,n) + t_cog  

                output_map[i].append(pulse)
                output_10[i].append(pulse2)

            radius[i] = np.array(vec)
            radius_old[i] = np.array(vec_old)

        frame['All_radius'] = radius
        frame['All_radius_old'] = radius_old
        frame['All_pulses'] = output_map
        frame['All_10'] = output_10
        self.PushFrame(frame)

tray = I3Tray()

########## SERVICES FOR GULLIVER ##########

#------------------- LET'S RUN SOME MODULES!  ------------------

#**************************************************
#                 Reader and whatnot
#**************************************************

tray.AddService("I3GulliverMinuitFactory","Minuit")(
    ("MinuitPrintLevel",-2),
    ("FlatnessCheck",True),
    ("Algorithm","MIGRAD"),
    ("MaxIterations",50000),
    ("MinuitStrategy",2),
    ("Tolerance",0.1),
    )

tray.AddService("I3CurvatureSeedServiceFactory","CurvSeed")(
    ("SeedTrackName", "Laputop"), # This is also the default                                                                        
    ("A", 6e-4),            # This comes from the original IT-26 gausspar function                                                  
    ("N",9.9832),
    ("D",63.5775)
    )

tray.AddService("I3CurvatureParametrizationServiceFactory","CurvParam")(
    ("FreeA", True),
    ("MinA", 0.0),
    ("MaxA", 2e-3),
    ("StepsizeA", 1e-5)
    )

tray.AddService("I3CurvatureParametrizationServiceFactory","CurvParam2")(
    ("FreeN",True),
    ("MinN",0),
    ("MaxN",200.0),
    ("StepsizeN",2.0)
    )

tray.AddService("I3CurvatureParametrizationServiceFactory","CurvParam3")(
    ("FreeD",True),
    ("MinD",0),
    ("MaxD",500.0),
    ("StepsizeD",2.0)
    )

datareadoutName = 'All_pulses'
badtanksName= "BadDomsList"

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

tray.AddService("I3LaputopLikelihoodServiceFactory","ToprecLike2")(
    ("datareadout", datareadoutName),
    ("badtanks", badtanksName),
    ("ldf", ""),      # do NOT do the LDF (charge) likelihood                                                                          
    ("curvature","gaussparfree")      # yes, do the CURVATURE likelihood                                                               
    )


tray.AddModule("I3LaputopFitter","CurvatureOnly")(
    ("SeedService","CurvSeed"),
    ("NSteps",3),                    # <--- tells it how many services to look for and perform                                          
    ("Parametrization1","CurvParam"),
    ("Parametrization2","CurvParam2"),
    ("Parametrization3","CurvParam3"),
    ("StoragePolicy","OnlyBestFit"),
    ("Minimizer","Minuit"),
    ("LogLikelihoodService","ToprecLike2"),     # the three likelihoods                                                                 
    ("LDFFunctions",["","",""]),   # do NOT do the LDF (charge) likelihood                                                              
    ("CurvFunctions",["gaussparfree","gaussparfree","gaussparfree"]) # yes, do the CURVATURE likelihood                                 
    )

tray.AddModule("I3Writer","EventWriter")(
    ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
    ("Filename",I3_OUTFILE),
)

wanted_inice_reco=["Millipede",
                           "MillipedeFitParams",
                           "Millipede_dEdX",
                           "Stoch_Reco",
                           "Stoch_Reco2",
                           "I3MuonEnergyLaputopCascadeParams",
                           "I3MuonEnergyLaputopParams"
                           ]

wanted_inice_cuts=['IT73AnalysisInIceQualityCuts']

wanted_general = ['I3EventHeader',
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
                  'Laputop',
                  'LaputopParams',
                  'All_pulses',
                  'All_radius',
                  'All_radius_old',
                  'ShowerCOG',
                  'CurvatureOnly',
                  'CurvatureOnlyParams'
              ]


## Output root file
root = I3ROOTTableService(ROOTFILE,"aTree")
tray.AddModule(I3TableWriter,'writer')(
    ("tableservice", root),
    ("keys",wanted_general+wanted_inice_reco+wanted_inice_cuts),
    ("subeventstreams",['InIceSplit',"ice_top"])
)



   
# Execute the Tray
# Just to make sure it's working!
tray.Execute()

