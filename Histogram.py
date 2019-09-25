import uproot as up
import numpy as np
import matplotlib.pylab as plt
import argparse

parser = argparse.ArgumentParser(description='Process DNN')

parser.add_argument('-i',
                    dest = 'input_file',
                    help='This is the input root file')

parser.add_arguments('-o',
                     dest = 'output_base',
                     help='Base of output files for histograms')

args = parser.parse_args()
