#!/usr/bin/python

## This tray segment will perform a three-step Laputop fit, using the parameters used in the 3 year analyses (IC79-IC86.2012).

from icecube import icetray, gulliver, gulliver_modules, lilliput

@icetray.traysegment
def LaputopStandard(tray, name, 
                    pulses='CleanedHLCTankPulses',
                    excluded='ClusterCleaningExcludedStations',
                    ShowerCOGSeed='ShowerCOG',
                    ShowerPlaneSeed='ShowerPlane',
                    snowfactor=2.1,
                    If = lambda f: True):

    ## Some more defaults
    fixcore = False     # do NOT keep core fixed

    ########## SERVICES FOR GULLIVER ##########

    ## The "simple lambda" snowservice
    tray.AddService("I3SimpleSnowCorrectionServiceFactory",name+"SimpleSnow")(
        ("Lambda", snowfactor)
    )

    ## This one is the standard one.
    tray.AddService("I3GulliverMinuitFactory",name+"Minuit")(
        ("MinuitPrintLevel",-2),  
        ("FlatnessCheck",True),  
        ("Algorithm","SIMPLEX"),  
        ("MaxIterations",5000),
        ("MinuitStrategy",2),
        ("Tolerance",0.01),    
        )
    
    ## The Seed service
    tray.AddService("I3LaputopSeedServiceFactory",name+"ToprecSeed")(
        ("InCore", ShowerCOGSeed),
        ("A",4.823e-4),
        ("N",0),
        ("InPlane", ShowerPlaneSeed),
        ("Beta",2.6),                    # first guess for Beta
        ("InputPulses",pulses)  # this'll let it first-guess at S125 automatically
        )
    
    ## Step 1:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam2")(
        ("FixCore", fixcore),        
        ("FixTrackDir", True),
        ("IsBeta", True),
        ("MinBeta", 2.9),   ## From toprec... 2nd iteration (DLP, using beta)
        ("MaxBeta", 3.1),
        ("maxLogS125",8.0),        # Default is 6., be a bit safer, although should never happen to be this large
        ("VertexStepsize",10.0),   # The COG is very good for contained events, don't step too far
        ("LimitCoreBoxSize", 200.0) 
    )

    ## Step 2:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam3")(
        ("FixCore", fixcore),        
        ("FixTrackDir", False),      # FREE THE DIRECTION!
        ("FreeA",True),
        ("FreeN",False),
        ("FreeD",False),
        ("IsBeta", True),
        ("MinBeta", 2.0),   ## From toprec... 3rd iteration (DLP, using beta)
        ("MaxBeta", 4.0),
        ("LimitCoreBoxSize", 15.0),
        ("maxLogS125",8.0),                   
        ## Use these smaller stepsizes instead of the defaults:
        ("VertexStepsize",2.0),      # default is 20
        ("SStepsize", 0.045),        # default is 1
        ("BetaStepsize",0.15)        # default is 0.6    
        )
        
    ## Step 3:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam4")(
        ("FixCore", fixcore),        
        ("FixTrackDir", True),
        ("IsBeta", True),
        ("MinBeta", 0.0),   
        ("MaxBeta", 10.0),
        ("LimitCoreBoxSize", 45.0),
        ("maxLogS125",8.0),        
        ## Use these smaller stepsizes instead of the defaults:
        ("VertexStepsize", 1.0),     # default is 20
        ("SStepsize", 0.045),        # default is 1
        ("BetaStepsize",0.15)        # default is 0.6 
        )
    
    tray.AddService("I3LaputopLikelihoodServiceFactory",name+"ToprecLike2")(
        ("DataReadout", pulses),
        ("BadTanks", excluded),
        ("DynamicCoreTreatment", 5.0),     # do the 5-meter core cut (used 11m before, but TF 5m and didn't see a big difference in the plots)
        ("SaturationLikelihood", True),
        ("MaxIntraStationTimeDiff",80.0),    # Don't use time fluctuating tanks for timing fits, could really mess up the hard work
        ("Curvature","gaussparfree"),      # NO timing likelihood (at first; this will be overridden)
        ("SnowServiceName",name+"SimpleSnow"),
        ("OldXYZ", True)  # For backward-compatibility for L3: DOM coordinates
        )
        
    ################# GULLIVERIZED FITTER MODULE #######################
    
    ## This module performs the three steps
    tray.AddModule("I3LaputopFitter",name)(
        ("SeedService",name+"ToprecSeed"),
        ("NSteps",3),            # <--- tells it how many services to look for and perform
        ("Parametrization1",name+"ToprecParam2"),   # the three parametrizations
        ("Parametrization2",name+"ToprecParam3"),
        ("Parametrization3",name+"ToprecParam4"),
        ("StoragePolicy","OnlyBestFit"),
        ("Minimizer",name+"Minuit"),
        ("LogLikelihoodService",name+"ToprecLike2"),     # the three likelihoods
        ("LDFFunctions",["dlp","dlp","dlp"]),
        ("CurvFunctions",["","gausspar","gausspar"]),   # VERY IMPORTANT : use time Llh for step 3, but fix direction!
        ("If",If)
        )
    

