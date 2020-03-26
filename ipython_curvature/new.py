#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:33:10 2020

@author: andy
"""

from icecube import *
from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, dataio, WaveCalibrator, wavedeform, wavereform
from icecube.payload_parsing import I3DOMLaunchExtractor
from icecube.icetray import I3Module
from icecube.icetop_Level3_scripts.segments import IceTopQualityCuts
import numpy as np
from icecube.dataclasses import *

infile = ['Subset_0.i3']
gain_doms = {}

class PutInts(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)

    def DAQ(self,frame):
        status = frame['I3DetectorStatus'].dom_status
        for i in frame['I3DetectorStatus'].dom_status.keys():
            gain_doms[i] = str(status[i].dom_gain_type)

    def Physics(self, frame):
        MCPE = frame['I3MCTree']
        for i in list(MCPE):
            print(i.ID)
        count = 0
        #for i in MCPE.keys():
        #    if count == 0:
        #        print(MCPE[i])
        #    count+=1
        
        self.PushFrame(frame)                 # push the frame
        


tray = I3Tray()


tray.AddModule("I3Reader","reader")(
        ("FileNameList", infile))

tray.AddModule(PutInts)

'''tray.AddModule("I3Writer","EventWriter")(
        ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
        ("Filename",'output_waveforms.i3.bz2'),
)
'''
tray.Execute(5)

np.save('DOM_gain',gain_doms)
