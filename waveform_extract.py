#!/home/amedina/build_stable/bin/python
import argparse
import os
import sys,getopt
#import logging
from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, frame_object_diff,simclasses,WaveCalibrator
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
HOME = '/data/user/amedina/Waveform_file/'

data_set_number = args.directory_number
directory = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/'%(data_set_number)
directory_list = os.listdir('/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/'%(data_set_number))
file_list = []
for i in directory_list:
    x = np.array(os.listdir(directory+i))
    y = [directory+i+'/'+m for m in x]
    file_list.append(y)

file_list_all = np.hstack(file_list)

file_list_all = np.array(np.array_split(file_list_all,200))
    
#print(file_list_all)
#### PUT YOUR FAVORITE GCD AND INPUT FILE HERE
# This particular example file lives in Madison.
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


def output_i3_root(i):
    infile = file_list_all[i]
    infile = [str(j) for j in infile]
    print(infile,i)

    #### NAME YOUR OUTPUT FILES HERE
    # Put this someplace that exists!  (This particular one is for KR)
    I3_OUTFILE = HOME+'/'+data_set_number+'/'+'Subset_%s'%(i)+".i3"
    ROOTFILE = HOME+'/'+data_set_number+'/'+'Subsets_%s'%(i)+".root"


    tray = I3Tray()

    ########## SERVICES FOR GULLIVER ##########

    datareadoutName="IceTopHLCSeedRTPulses"
    badtanksName= "IceTopHLCSeedRTExcludedTanks"

    ## The "simple lambda" snowservice
    #tray.AddService("I3SimpleSnowCorrectionServiceFactory","SimpleSnow21")(
   #    ("Lambda", 2.1)
    #    )

    ## This one is the standard one.

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


    tray.AddService("I3LaputopLikelihoodServiceFactory","ToprecLike2")(
        ("datareadout", datareadoutName),
        ("badtanks", badtanksName),
        ("ldf", ""),      # do NOT do the LDF (charge) likelihood
        ("curvature","gaussparfree")      # yes, do the CURVATURE likelihood
        #    ("SnowServiceName","SimpleSnow21")
    )


    #------------------- LET'S RUN SOME MODULES!  ------------------

    #**************************************************
    #                 Reader and whatnot
    #**************************************************
    tray.AddModule("I3Reader","reader")(
        ("FileNameList", [GCDfile]+infile)
    )

    tray.AddModule(Remove_WaveformRange)

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

    from icecube.icetop_Level3_scripts.functions import count_stations
    from icecube.icetop_Level3_scripts.modules import FilterWaveforms
    from icecube.icetop_Level3_scripts.segments.extract_waveforms import ExtractWaveforms

    tray.AddModule(FilterWaveforms, 'FilterWaveforms',   #Puts IceTopWaveformWeight in the frame.                                                     
                   pulses=icetop_globals.icetop_hlc_pulses,
                   If = lambda frame: icetop_globals.icetop_hlc_pulses in frame and count_stations(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, icetop_globals.icetop_hlc_pulses)) >= 5)
    tray.AddSegment(ExtractWaveforms, 'IceTop')
                    #If= lambda frame: "IceTopWaveformWeight" in frame and frame["IceTopWaveformWeight"].value!=0)


    itpulses='IceTopHLCSeedRTPulses'
    
    tray.AddSegment(level3_IceTop, "level3_IceTop",
                    detector='IC86.2012',
                    do_select = False,
                    isMC=True,
                    add_jitter=False,
                    snowLambda=None
                    )

    def fix_spe(frame,pulses):
        if pulses in frame: # Although InIcePulses should always be there                                                                                                      

            corr = dataclasses.I3RecoPulseSeriesMapApplySPECorrection(pulses,"I3Calibration")
            corrpulses = corr.apply(frame)
            frame.Delete(pulses)
            frame[pulses] = corrpulses
            return True

    tray.AddModule(fix_spe, "fixspe",pulses=icetop_globals.names["Level3"]["InIcePulses"],Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

    tray.AddSegment(level3_Coinc,"level3_Coinc",
                    Detector='IC86.2012',
                    isMC=True,
                    do_select=False,
                    IceTopTrack='Laputop',
                    IceTopPulses=itpulses,
    )

    #Keys to Keep

    wanted_general=['I3EventHeader',
                    icetop_globals.filtermask,
                    'I3TriggerHierarchy']


    wanted_general+=['MCPrimary',
                     'MCPrimaryInfo',
                     'AirShowerComponents',
                     'IceTopComponentPulses_Electron',
                     'IceTopComponentPulses_ElectronFromChargedMesons',
                     'IceTopComponentPulses_Gamma',
                     'IceTopComponentPulses_GammaFromChargedMesons',
                     'IceTopComponentPulses_Muon',
                     'IceTopComponentPulses_Hadron',
    ]

    wanted_icetop_filter=['IceTop_EventPrescale',
                          'IceTop_StandardFilter',
                          'IceTop_InFillFilter']
 
    wanted_icetop_pulses=[icetop_globals.icetop_hlc_pulses,
                          icetop_globals.icetop_slc_pulses,
                          icetop_globals.icetop_clean_hlc_pulses,
                          icetop_globals.icetop_tank_pulse_merger_excluded_tanks,
                          icetop_globals.icetop_cluster_cleaning_excluded_tanks,
                          icetop_globals.icetop_HLCseed_clean_hlc_pulses,
                          icetop_globals.icetop_HLCseed_excluded_tanks,
                          icetop_globals.icetop_HLCseed_clean_hlc_pulses+'_SnowCorrected',
                          'TankPulseMergerExcludedSLCTanks',
                          'IceTopLaputopSeededSelectedHLC',  
                          'IceTopLaputopSeededSelectedSLC',                                                                                                                                               
                          'IceTopLaputopSmallSeededSelectedHLC',  
                          'IceTopLaputopSmallSeededSelectedSLC',                                                                                                                                         
                          ]

    wanted_icetop_waveforms=['IceTopVEMCalibratedWaveforms',
                             'IceTopWaveformWeight',
                             'ReextractedHLCFADCWaveforms',
                             'ReextractedSLCWaveforms']

    wanted_icetop_reco=['ShowerCOG',
                        'ShowerPlane',
                        'ShowerPlaneParams',
                        'Laputop',
                        'LaputopParams',
                        'LaputopSnowDiagnostics',
                        'LaputopSmall',
                        'LaputopSmallParams',
                        'IsSmallShower'
       #                 "CurvatureOnly",
       #                 "CurvatureOnlyParams"
                        ]
    
    wanted_icetop_cuts=['Laputop_FractionContainment',
                        'Laputop_OnionContainment',
                        'Laputop_NearestStationIsInfill',
                        'StationDensity',
                        'IceTopMaxSignal',
                        'IceTopMaxSignalInEdge',
                        'IceTopMaxSignalTank',
                        'IceTopMaxSignalString',
                        'IceTopNeighbourMaxSignal',
                        'IT73AnalysisIceTopQualityCuts',
                        ]

    wanted=wanted_general+wanted_icetop_filter+wanted_icetop_pulses+wanted_icetop_waveforms+wanted_icetop_reco+wanted_icetop_cuts
    
    wanted_sta2 = wanted_general + wanted_icetop_filter + \
                    [icetop_globals.icetop_hlc_pulses,
                    icetop_globals.icetop_slc_pulses,
                    icetop_globals.icetop_clean_hlc_pulses,
                    icetop_globals.icetop_cluster_cleaning_excluded_tanks,
                    icetop_globals.icetop_tank_pulse_merger_excluded_tanks,
                    'TankPulseMergerExcludedSLCTanks']

    #InIce Keys to Keep
    wanted_inice_pulses=[icetop_globals.inice_pulses,
                         icetop_globals.inice_coinc_pulses,
                         icetop_globals.inice_clean_coinc_pulses,
                         icetop_globals.inice_clean_coinc_pulses+"TimeRange",
                         icetop_globals.inice_clean_coinc_pulses+"_Balloon",
                         "SaturationWindows",
                         "CalibrationErrata",
                         'SRT'+icetop_globals.inice_coinc_pulses,
                         'NCh_'+icetop_globals.inice_clean_coinc_pulses]

    wanted_inice_reco=["Millipede",
                       "MillipedeFitParams",
                       "Millipede_dEdX",
                       "Stoch_Reco",
                       "Stoch_Reco2",
                       "I3MuonEnergyLaputopCascadeParams",
                       "I3MuonEnergyLaputopParams"
    ]
   
    wanted_inice_cuts=['IT73AnalysisInIceQualityCuts']

    wanted_inice_muon=['CoincMuonReco_LineFit',
                       'CoincMuonReco_SPEFit2',
                       'CoincMuonReco_LineFitParams',
                       'CoincMuonReco_SPEFit2FitParams',
                       'CoincMuonReco_MPEFit',
                       'CoincMuonReco_MPEFitFitParams',
                       'CoincMuonReco_MPEFitMuEX',
                       'CoincMuonReco_CVMultiplicity',
                       'CoincMuonReco_CVStatistics',
                       'CoincMuonReco_MPEFitCharacteristics',
                       'CoincMuonReco_SPEFit2Characteristics',
                       'CoincMuonReco_MPEFitTruncated_BINS_Muon',
                       'CoincMuonReco_MPEFitTruncated_AllBINS_Muon',
                       'CoincMuonReco_MPEFitTruncated_ORIG_Muon',
                       'CoincMuonReco_SPEFit2_D4R_CascadeParams',
                       'CoincMuonReco_SPEFit2_D4R_Params',
                       'CoincMuonReco_MPEFitDirectHitsC'
    ]

    wanted=wanted+wanted_inice_pulses+wanted_inice_reco+wanted_inice_cuts+ wanted_inice_muon

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
        ("keys",[ "Laputop", 
                  "LaputopParams",
                  #"CurvatureOnly",
                  #"CurvatureOnlyParams",
                  "Millipede",
                  "MillipedeFitParams",
                  "Millipede_dEdX",
                  "Stoch_Reco",
                  "Stoch_Reco2",
                  "I3EventHeader",
                  "MCPrimary",
                  "MCPrimaryInfo",
                  "IT73AnalysisIceTopQualityCuts",
                  "IT73AnalysisInIceQualityCuts",
                  'IceTopVEMCalibratedWaveforms']),
        ("subeventstreams",['InIceSplit',"ice_top"])
    )

   
    # Execute the Tray
    # Just to make sure it's working!
    tray.Execute()

pool = mp.Pool(processes=5)
pool.map(output_i3_root,x_list)
