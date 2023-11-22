# CosmicRay Analysis

**Author:** Andres Medina

**Github:** [andyPhysics](https://github.com/andyPhysics)

**Institution:** The Ohio State University

## Overview

The CosmicRay analysis code, authored by Andres Medina, aims to investigate the composition of cosmic rays using simulation data. The code processes data located in the directory `/data/user/amedina/CosmicRay`. Below is an organized overview of the primary directories and files within the project.

## Directory Structure

### Simulation Data Analysis

1. **Analysis**: I3 files and root files with curvature and thickness fit.
2. **Analysis_data**: I3 files and root files with curvature and thickness fit for experimental data.
3. **csv_files**: CSV files with essential information for training networks.
4. **csv_files_data**: CSV files with necessary information for experimental data.
5. **I3_files**: I3 files with coincident events used for curvature and thickness calculations, including waveforms.
6. **I3_files_data**: I3 files for experimental data with coincident events, waveforms, and a 10% cut.

### Top-Level Files

#### Jupyter Notebooks

- **Check-fit.ipynb**: Used for checking the fit of thickness parameters.
- **Untitled.ipynb**: Extracts waveforms for verification.

#### Python Files

- **analysis.py**: Conducts curvature and thickness reconstruction.
- **analysis_data.py**: Performs curvature and thickness reconstruction for experimental data.
- **curvature.py**: Previous version of the curvature file.
- **get_event.py**: Extracts about 1000 events to produce Events.csv for fit checking.
- **waveform_extract_2.py**: Extracts necessary events with coincident events and waveforms.
- **waveform_extract_data.py**: Extracts necessary events with coincident events and waveforms for experimental data.
- **my_laputop_segment.py**: Laputop segment used for curvature reconstruction in analysis.
- **methods.py, methods2.py, vectorDiff.py**: Provide methods for changing basis and necessary thickness reconstruction variables.

#### Submitting Files

- **analysis.sh**: Runs the analysis.
- **waveform_extract.sh**: Obtains waveform files from I3Files.
- **waveform_extract_data.sh**: Obtains waveform files from I3Files_data.
- **output_files.sh**: Creates a list of files for analysis.

## Unused Section

### ipython_curvature

This section can be ignored.

## Neural Network (NN) Section

### Files

- **data.csv, Proton.csv, Helium.csv, Oxygen.csv, Iron.csv**: Key files for training and storing network models.

#### Jupyter Notebooks

- **Composition.ipynb**: Final notebook producing KDEs for composition comparisons.
- **model_plotting.ipynb**: Generates plots for different models (H3A, H4A, GST, GST-3gen).
- **NN_fitter.ipynb**: Tests individual networks for overfitting and reconstruction capabilities.

#### Python Files

- **extract_data.py**: Extracts parameters from the Analysis folder into CSV files in the csv_files folder.
- **extract_data2.py**: Similar to extract_data.py but for experimental data.
- **methods.py**: Provides methods used in NN and NN_fitter2, defining cuts, also used in Jupyter notebooks.
- **NN.py**: Old version used for training the network.
- **NN_fitter2.py**: Updated version, the model used to produce networks for reconstructing log10(E/GeV) and Xmax.

#### Submitting Files

- **submit.sh**: Used to extract the CSV files.
- **run_random_state.sh**: Runs a different seed for NN_fitter2.py to produce five distinct neural networks trained on different data sets.


