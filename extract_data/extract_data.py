#!/usr/bin/env python

import numpy as np
import h5py
import argparse
import uproot
from collections import OrderedDict
import itertools
import random
import datetime
import sys,os,fnmatch
from scipy.optimize import curve_fit
from scipy.stats import chisquare


