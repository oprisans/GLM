##Project Title
Generalized Linear Model Simulations


##Description
The data were recorded form mice both during exploration and REM sleep and included in the Excel file. The number of neurons is variable for each mouse and depends on experimental conditions (see C. Blanco-Centurion, S. Luo, D. J. Spergel, A. Vidal-Ortiz, S. A. Oprisan, A. N. Van den Pol, M. Liu, P. J. Shiromani, Dynamic network activation of hypothalamic mch neurons in rem sleep and exploratory behavior, J Neurosci 39 (25) (2019) 4986–4998. doi:10.1523/JNEUROSCI.0305-19.2019). The Matlab code reads the data and estimates the coupling coefficients using #y_#x GLM where y = predicted variable and x = predictor. The Matlab codes are for 1y_1x, hich is a simple correlation measure and 1y_Allx that produces one prection based on all th other recordings. Code for articial data generation is provided as a ground truth to check the pipeline. Use provided z-score calcium fluorescence data to estimate GLM couplings.

## How to cite
If you use this code/data, please cite:
Manuscript citation / preprint:
Zenodo release: 10.5281/zenodo.18102531

## Contents
a) Reproducible pipeline to generate: Figure/Table
b) Implementation in MATLAB
c) Includes: minimal example dataset
d) Does NOT include: raw miniscope movies / proprietary files

## Repository structure
raw_data_analysis - use it on raw z-score data provided
ground_truth - generates artificaial data from a known linear model thatcan then be run through the same piplene as the raw data
publication_scripts - agregaes data and generates statistics and figures

## Getting Started
##Installing:
1. Install MATLAB 2025b with ALL toolboxes
2. Clone this repo
3. Open MATLAB in repo root
4. Run individual files to analyze data or generate groubd truth using the same pipeline.
5. Outputs are written to: `results/` and `figures/`

## Authors: Dr. Srinel A. Oprisan
oprisans@cofc.edu

## License:
This project is licensed under the GLP License - see the LICENSE.md file for details
Data license (often CC-BY 4.0)
