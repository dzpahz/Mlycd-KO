# README

## Description
This folder contains three R scripts used for the analysis of proteomics and metabolomics data described in the paper "MCD Enzyme Deficiency Drives Early Mortality and Multiorgan Metabolic Disruption in Dahl Salt-Sensitive Rats".

## Files
1. **`DEP process.R`**  
   Performs data processing, normalization, and statistical analysis of proteomics datasets exported from Spectronaut.  
   Main outputs: differential protein lists.  

2. **`metabolomics.R`**  
   Handles raw metabolite data processing from Skyline transitions.  
   Non-lysine and non-TCA metabolites were excluded during data processing.  

3. **`metabolite analysis and heatmap.R`**  
   Conducts metabolomics data analysis and visualization.  
   Generates heatmaps.  

## Dependencies
- Required packages: `dplyr`, `rstatix`, `ggplot2`, `tidyr`, etc.  

## Note
- Proteomics data have been deposited in the ProteomeXchange Consortium under accession number **PXD062350**.
- Metabolomics data are available in the **MassIVE** repository; details can be found in the paper.
