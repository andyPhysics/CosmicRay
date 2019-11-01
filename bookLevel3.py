'''
Example booking script for level3 IT files.
Not all keys are listed here, you should add your own selection...
'''
import glob, re
from argparse import ArgumentParser
from icecube import icetray, dataio, dataclasses, phys_services
from icecube.tableio import I3TableWriter
from icecube.rootwriter import I3ROOTTableService
from icecube.icetop_Level3_scripts import icetop_globals
from icecube.frame_object_diff.segments import uncompress
from I3Tray import *


parser = ArgumentParser(usage='%s [options] -o <filename>.i3[.bz2|.gz] {i3 file list}'%os.path.basename(sys.argv[0]))
parser.add_argument("-o", "--output", action="store", type=str, dest="output", help="Output file name", metavar="BASENAME")
parser.add_argument('-i','--input',action='store',dest='inputFiles',help="Input file(s)")

(args) = parser.parse_args()

# Check whether output is okay.
if not args.output:
    icetray.logging.log_error("Output file not specified!")
    ok=False
else:
    if args.output[-5:]==".root":
        table_service = I3ROOTTableService(args.output)
    elif args.output[-3:]==".h5":
        icetray.logging.log_error("I do not know how to handle h5 files yet.")
        ok=False
    else:
        icetray.logging.log_error("Wrong extension for booking.")
        ok=False


tray=I3Tray()
tray.AddModule("I3Reader","reader", FilenameList=[args.inputFiles])

it_gen=['I3EventHeader']

it_gen=it_gen+['MCPrimary','MCPrimaryInfo']

it_filter=['IceTop_EventPrescale',
           'IceTop_StandardFilter',
           'IceTop_InFillFilter']

it_pulses=[dict(key =icetop_globals.icetop_clean_hlc_pulses,
                converter = dataclasses.converters.I3RecoPulseSeriesMapConverter(bookGeometry=True),
                name = icetop_globals.icetop_clean_hlc_pulses),
           dict(key =icetop_globals.icetop_clean_hlc_pulses,
                converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(),
                name = icetop_globals.icetop_clean_hlc_pulses+'_info'),
           dict(key =icetop_globals.icetop_HLCseed_clean_hlc_pulses,
                converter = dataclasses.converters.I3RecoPulseSeriesMapConverter(bookGeometry=True),
                name = icetop_globals.icetop_HLCseed_clean_hlc_pulses),
           dict(key =icetop_globals.icetop_HLCseed_clean_hlc_pulses,
                converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(),
                name = icetop_globals.icetop_HLCseed_clean_hlc_pulses+'_info'),
           dict(key =icetop_globals.icetop_HLCseed_clean_hlc_pulses+"_SnowCorrected",
                converter = dataclasses.converters.I3RecoPulseSeriesMapConverter(bookGeometry=True),
                name = icetop_globals.icetop_HLCseed_clean_hlc_pulses+"_SnowCorrected"),
           dict(key =icetop_globals.icetop_HLCseed_clean_hlc_pulses+"_SnowCorrected",
                converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(),
                name = icetop_globals.icetop_HLCseed_clean_hlc_pulses+'_SnowCorrected_info')
           ]


it_reco=["Laputop","LaputopParams",
         dict(key='Laputop',
              converter = phys_services.converters.I3RecoInfoConverter('IceTopHLCSeedRTPulses'),
              name = 'Laputop_info'),
         ]

it_cuts=["IT73AnalysisIceTopQualityCuts"]


book_keys=it_gen+it_filter+it_pulses+it_reco+it_cuts
from icecube import ddddr, common_variables
# Make some selection
ic_pulses=[dict(key ='SRT'+icetop_globals.inice_coinc_pulses,
                converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(),
                name = 'SRT'+icetop_globals.inice_coinc_pulses+'_info'),
           dict(key =icetop_globals.inice_clean_coinc_pulses,
                converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(),
                name = icetop_globals.inice_clean_coinc_pulses+'_info')]
   

# Also some selection
ic_reco=["Millipede",
         "MillipedeFitParams",
         "Millipede_dEdX",
         "I3MuonEnergyLaputopParams",
         'CoincMuonReco_MPEFit',
         'CoincMuonReco_MPEFitFitParams',
         'CoincMuonReco_SPEFit2_D4R_Params',
         "Stoch_Reco",
         "Stoch_Reco2",
         "CurvatureOnly",
         "CurvatureOnlyParams"]
    
ic_cuts=['IT73AnalysisInIceQualityCuts']

book_keys=book_keys+ic_pulses+ic_reco+ic_cuts
    
tray.AddModule(I3TableWriter, "table_writer",
               TableService = table_service,
               SubEventStreams = ['ice_top'],
               Keys = book_keys)

# Execute the Tray                                                                                                                                                                                       

tray.Execute()
    

