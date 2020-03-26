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

infile = ['/home/andy/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz']

geometry = {}

class PutInts(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("Geometry",
                          "Positions of DOMs")

    def Configure(self):
        self.geometry = self.GetParameter('Geometry')

    def Geometry(self,frame):
        for keys in frame["I3Geometry"].omgeo.keys():
            self.geometry[keys] = frame["I3Geometry"].omgeo[keys].position
        self.PushFrame(frame)



tray = I3Tray()


tray.AddModule("I3Reader","reader")(
        ("FileNameList", infile))

tray.AddModule(PutInts,Geometry = geometry)

'''tray.AddModule("I3Writer","EventWriter")(
        ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
        ("Filename",'output_waveforms.i3.bz2'),
)
'''
tray.Execute()

np.save('geometry',{'geometry':geometry})

