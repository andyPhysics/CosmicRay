# CosmicRay Analysis

**Author:** Andres Medina

**Github:** [andyPhysics](https://github.com/andyPhysics)

**Institution:** The Ohio State University

## Overview

Welcome to the CosmicRay analysis code, a project by Andres Medina. This code aims to investigate the composition of cosmic rays using simulation data. The project is organized into various directories, each serving a specific purpose and containing essential files. The data is located in the directory `/data/user/amedina/CosmicRay`. Below is a structured overview of the primary directories and files within the project.

## Directory Structure

### Simulation Data Analysis

1. **NN**: Neural Network files and scripts for composition analysis.
2. **analysis_scripts**: Python scripts for curvature and thickness reconstruction.
3. **data**: CSV files with essential information for training networks.
4. **data_files**: CSV files with necessary information for experimental data.
5. **I3_files**: I3 files with coincident events used for curvature and thickness calculations, including waveforms.
6. **I3_files_data**: I3 files for experimental data with coincident events, waveforms, and a 10% cut.

### Top-Level Files

#### Jupyter Notebooks

- **notebooks/Check-fit.ipynb**: Utilized for checking the fit of thickness parameters.
- **notebooks/Untitled.ipynb**: Extracts waveforms for verification.

#### Python Files

- **scripts/analysis.py**: Conducts curvature and thickness reconstruction.
- **scripts/analysis_data.py**: Performs curvature and thickness reconstruction for experimental data.
- **scripts/curvature.py**: Previous version of the curvature file.
- **scripts/get_event.py**: Extracts about 1000 events to produce Events.csv for fit checking.
- **scripts/waveform_extract_2.py**: Extracts necessary events with coincident events and waveforms.
- **scripts/waveform_extract_data.py**: Extracts necessary events with coincident events and waveforms for experimental data.
- **scripts/my_laputop_segment.py**: Laputop segment used for curvature reconstruction in analysis.
- **scripts/methods.py, scripts/methods2.py, scripts/vectorDiff.py**: Provide methods for changing basis and necessary thickness reconstruction variables.

#### Submitting Files

- **analysis.sh**: Runs the analysis.
- **scripts/waveform_extract.sh**: Obtains waveform files from I3Files.
- **scripts/waveform_extract_data.sh**: Obtains waveform files from I3Files_data.
- **scripts/output_files.sh**: Creates a list of files for analysis.

## Neural Network (NN) Section

### Files

- **NN/data.csv, NN/Proton.csv, NN/Helium.csv, NN/Oxygen.csv, NN/Iron.csv**: Key files for training and storing network models.

#### Jupyter Notebooks

- **NN/notebooks/Composition.ipynb**: Final notebook producing KDEs for composition comparisons.
- **NN/notebooks/model_plotting.ipynb**: Generates plots for different models (H3A, H4A, GST, GST-3gen).
- **NN/notebooks/NN_fitter.ipynb**: Tests individual networks for overfitting and reconstruction capabilities.

#### Python Files

- **NN/scripts/extract_data.py**: Extracts parameters from the analysis_scripts folder into CSV files in the data folder.
- **NN/scripts/extract_data2.py**: Similar to extract_data.py but for experimental data.
- **NN/scripts/methods.py**: Provides methods used in NN and NN_fitter2, defining cuts, also used in Jupyter notebooks.
- **NN/scripts/NN_fitter2.py**: Updated version, the model used to produce networks for reconstructing log10(E/GeV) and Xmax.

#### Submitting Files

- **NN/submit.sh**: Used to extract the CSV files.
- **NN/run_random_state.sh**: Runs a different seed for NN_fitter2.py to produce five distinct neural networks trained on different data sets.


