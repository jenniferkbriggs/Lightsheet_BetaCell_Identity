# README

This repository contains code for the analysis of 3D lightsheet calcium recordings of pancreatic islets. The findings from this analysis are published in:

**Jin Erli, Briggs Jennifer K., Benninger Richard K.P., Merrins Matthew J. (2024)**  
*Glucokinase activity controls subpopulations of β-cells that alternately lead islet Ca2+ oscillations*  
eLife 13:RP103068  
[https://doi.org/10.7554/eLife.103068.1](https://doi.org/10.7554/eLife.103068.1)

**Code:** Jennifer Briggs (2022)  
**Data:** Erli Jin (2021-2022)  
**Principal Investigators:** Richard Benninger, Matthew Merrins

## Repository Structure

The repository is organized into the following directories:

### 1. Run Scripts

- **RunNetworkandWave.m:** Main script for performing consistency analysis as presented in the publication.
- **RunNetworkandWave_3Dvs2D.m:** Similar to *RunNetworkandWave.m* but includes comparisons between 3D and 2D analyses.

### 2. AnalysisFunctions

- Contains all functions used for network and wave analysis, including the assessment of subpopulation consistency. These functions are utilized by the run scripts.

### 3. AnalyzeStructures

- Scripts in this folder process the output structures generated by the run scripts. 
- They parse the data and export results as CSV files for statistical analysis and plotting.

## Getting Started

1. **Prerequisites:**
   - Ensure MATLAB is installed with the necessary toolboxes for matrix computations and data visualization.
   - Clone and set up the following repository for additional functional and structural network analysis:  
     [https://github.com/jenniferkbriggs/Functional_and_Structural_Networks.git](https://github.com/jenniferkbriggs/Functional_and_Structural_Networks.git)
2. **Running Analyses:**
   - Execute *RunNetworkandWave.m* for standard consistency analysis.
   - Execute *RunNetworkandWave_3Dvs2D.m* for comparative analysis between 3D and 2D data.
3. **Post-Processing:** Use the scripts in the *AnalyzeStructures* folder to convert analysis results into CSV files for further statistical evaluation.

## Contact

For questions or issues, please contact Jennifer Briggs or the principal investigators, Richard Benninger and Matthew Merrins.

