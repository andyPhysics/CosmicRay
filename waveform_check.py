#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:33:10 2020

@author: andy
"""
from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, dataio, WaveCalibrator, wavedeform, wavereform
from icecube.payload_parsing import I3DOMLaunchExtractor
from icecube.icetray import I3Module
from icecube.icetop_Level3_scripts.segments import level3_IceTop, level3_Coinc, IceTopQualityCuts


GCD = 'GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz'
infile = [GCD,'Level2_IC86_corsika_icetop.010410.000333.i3.bz2']

tray = I3Tray()


class Remove_WaveformRange(icetray.I3Module):
    def DAQ(self,frame):
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


  

tray.AddModule("I3Reader","reader")(
        ("FileNameList", infile))

tray.AddModule(Remove_WaveformRange)

tray.AddModule("I3WaveCalibrator", "sedan",
    Launches="CleanIceTopRawData",
    Waveforms="CalibratedWaveform",
    Errata="BorkedOMs",
    ATWDSaturationMargin=123,
    FADCSaturationMargin=0,
    )

tray.AddModule("I3WaveformSplitter", "splitter",
    Input="CalibratedWaveform",
    HLC_ATWD="CalibratedATWD",
    HLC_FADC="CalibratedFADC_HLC",
    SLC="CalibratedFADC_SLC",
    PickUnsaturatedATWD=True,
    )

tray.AddSegment(level3_Coinc,
                 Detector='IC86',
                 IceTopTrack='Laputop',
                 InIcePulses='InIcePulses',
                 IceTopPulses='CleanedHLCTankPulses',
                 isMC=True,
                 do_select=False,
                 domeff=0.99, #Beware!                                                                                                                                                                
                 spline_dir="/data/sim/sim-new/downloads/spline-tables/"
                 )

tray.AddSegment(IceTopQualityCuts,
                      detector='IC86',
                      pulses='CleanedHLCTankPulses',
                      isMC=True,
                      reco_track='Laputop',
                      removeOrNot=False)



tray.AddModule("I3Writer","EventWriter")(
        ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
        ("Filename",'output_waveforms.i3.bz2'),
    )

tray.Execute()
