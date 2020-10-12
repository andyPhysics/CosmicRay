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

    tray.AddService("I3GulliverMinuitFactory",name+"Minuit")(
        ("MinuitPrintLevel",-2),
        ("FlatnessCheck",True),
        ("Algorithm","MIGRAD"),
        ("MaxIterations",50000),
        ("MinuitStrategy",2),
        ("Tolerance",0.1),
    )

    tray.AddService("I3CurvatureSeedServiceFactory",name+"CurvSeed")(
        ("SeedTrackName", "Laputop"), # This is also the default                                                                                                                     
        ("A", 6e-4),            # This comes from the original IT-26 gausspar function                                                                                               
        ("N",9.9832),
        ("D",63.5775)
    )

    tray.AddService("I3CurvatureParametrizationServiceFactory",name+"CurvParam")(
        ("FreeA", True),
        ("MinA", 0.0),
        ("MaxA", 2e-3),
        ("StepsizeA", 1e-5)
    )

    tray.AddService("I3CurvatureParametrizationServiceFactory",name+"CurvParam2")(
        ("FreeN",True),
        ("MinN",0),
        ("MaxN",200.0),
        ("StepsizeN",2.0)
    )

    tray.AddService("I3CurvatureParametrizationServiceFactory",name+"CurvParam3")(
        ("FreeD",True),
        ("MinD",0),
        ("MaxD",500.0),
        ("StepsizeD",2.0)
    )


    tray.AddService("I3LaputopLikelihoodServiceFactory",name+"ToprecLike2")(
        ("datareadout", pulses),
        ("badtanks", excluded),
        ("ldf", ""),      # do NOT do the LDF (charge) likelihood                                                                                                                    
        ("curvature","gaussparfree")      # yes, do the CURVATURE likelihood                                                                                                         
        #    ("SnowServiceName","SimpleSnow21")                                                                                                                                      
    )

    tray.AddModule("I3LaputopFitter","CurvatureOnly")(
        ("SeedService",name+"CurvSeed"),
        ("NSteps",1),                    # <--- tells it how many services to look for and perform                                                                                   
        ("Parametrization1",name+"CurvParam"),
        #("Parametrization2",name+"CurvParam2"),
        #("Parametrization3",name+"CurvParam3"),
        ("StoragePolicy","OnlyBestFit"),
        ("Minimizer",name+"Minuit"),
        ("LogLikelihoodService",name+"ToprecLike2"),     # the three likelihoods                                                                                                          
        ("LDFFunctions",["","",""]),   # do NOT do the LDF (charge) likelihood                                                                                                       
        ("CurvFunctions",["gaussparfree","gaussparfree","gaussparfree"]) # yes, do the CURVATURE likelihood                                                                          
    )
    

