# CosmicRay
Analysis of CosmicRay Composition
Author: Andres Medina
Github: andyPhysics
Institution: The Ohio State University

This code takes works with simulation data found in the following directory:
/data/user/amedina/CosmicRay

You will find the following directories:

Analysis -- I3files and root files with curvature and thickness fit

Analysis_data -- I3files and root files with curvature and thickness fit for experimental data

csv_files -- csv files containing all the necessary information for training networks

csv_files_data -- csv files containing all the necessary information for the experimental data

I3_files -- I3files that only contains coincident events, these are the files that are use to
 calculate the curvature and thickness. Contains also waveforms
 
I3_files_data -- I3files that only contains coincident events, these are the files that are used to calculate curvature and thickness for experimentald data. Contains waveforms and is a 10\% cut. 

## Top Level
 
Majority of the top level files are those used for data comparisons and extracting the necessary variables and reconstructions

###ipython notebooks

Check-fit.ipynb -- File used to check the fit for thickness parameters

Untitled.ipynb -- Extracts a waveform, for checks

### python files

Main files for the analysis to extract and reconstruct the files:

analysis.py -- Analysis file that performs the curvature and thickness reconstruction

analysis_data.py -- Analysis file that performs the curvature and thickness reconstruction for experimental data

curvature.py -- previous version of curvature

get_event.py -- file used to extract about 1000 events to produce Events.csv. This is for checking fits of the thickness parameters done in Check-fit.ipynb

waveform_extract_2.py -- File used to extract necessary events and only those that have coincident events and have waveforms. 

waveform_extract_data.py -- File used to extract necessary events and only those that have coincident events and have waveforms for experimental data

my_laputop_segment.py -- laputop segment used for curvature reconstruction in analysis

The three files below provide the methods for changing basis and those necessary for thickness reconstruction variables. 

methods.py 

methods2.py

vectorDiff.py

### Submitting files

analysis.sh -- file to run analysis

waveform_extract.sh -- obtain waveform files which are those in I3Files above

waveform_extract_data.sh -- obtain waveform files which are those in I3Files_data above

output_files.sh -- just creates list of files to be analyzed. 

Most of these files are run either in cobalt or you could run them by ./ > filename

## ipython_curvature

Feel free to ignore

## NN

**files** -- This contains all the files needed to train the network, along with the models that I have trained

Key files to keep in mind are data.csv, Proton.csv, Helium.csv, Oxygen.csv, and Iron.csv

Models - 

- NN -- model_coinc_best_[0-4].h5
- Decision Tree -- Energy_model.pkl
- Linear Model -- Xmax_bias_correction_[0-4].pkl

Correction number correspond to model_coinc_best number

### ipython Notebooks

Composition.ipynb -- This is the final notebook that takes networks and produces the KDEs for composition comparisons

model_plotting.ipynb -- Produces plots for the different models compared. This includes H3A, H4A, GST, and GST-3gen

NN_fitter.ipynb -- This file is used for testing individual networks to ensure that there is no overfitting and it has the capabilities of reconstructing both energy and Xmax.

## python files

extract_data.py -- This file extracts all of the necessary parameters from the Analysis folder into csv files that are stored in the csv_files folder. These files that are processed are stored in the files folder for quick access.

extract_data2.py -- This is similar to above but for data.

methods.py -- Provides the majority of the methods used in NN and NN_fitter2. This is where you can define the cuts. This is also used in the ipython notebooks. 

NN.py -- Old version used for training network

NN_fitter2.py -- Updated version and was the model used to produce the networks for reconstructing log10(E/GeV) and Xmax.

### Submitting files

submit.sh -- Used to extract the csv files. 

run_random_state.sh -- Runs a different seed for NN_fitter2.py to produce 5 distinct neural networks trained on a different set of the data. 



