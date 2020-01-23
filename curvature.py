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
HOME = '/data/user/amedina/CosmicRay/Curvature'

data_set_number = args.directory_number
directory = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/'
file_list = os.listdir(directory+data_set_number)
file_list_all = np.array([directory+data_set_number+'/'+i for i in file_list])
file_list_all = np.array(np.array_split(file_list_all,200))
#print(file_list_all)
#### PUT YOUR FAVORITE GCD AND INPUT FILE HERE
# This particular example file lives in Madison.
GCDfile = '/data/sim/IceTop/2012/filtered/CORSIKA-ice-top/%s/level2/0000000-0000999/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz'%(data_set_number)
x_list = list(range(file_list_all.shape[0]))
x_list = np.array(x_list)

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
 
    ## Output root file
    root = I3ROOTTableService(ROOTFILE,"aTree")
    tray.AddModule(I3TableWriter,'writer')(
        ("tableservice", root),
        ("keys",[ "Laputop", 
                  "LaputopParams",
                  "CurvatureOnly",
                  "CurvatureOnlyParams",
                  "Millipede",
                  "MillipedeFitParams",
                  "Millipede_dEdX",
                  "Stoch_Reco",
                  "Stoch_Reco2",
                  "I3EventHeader",
                  "MCPrimary",
                  "MCPrimaryInfo",
                  "IT73AnalysisIceTopQualityCuts",
                  "IT73AnalysisInIceQualityCuts"]),
        ("subeventstreams",["ice_top"])
    )
 
   
    # Execute the Tray
    # Just to make sure it's working!
    tray.Execute()

pool = mp.Pool(processes=10)
pool.map(output_i3_root,x_list)

