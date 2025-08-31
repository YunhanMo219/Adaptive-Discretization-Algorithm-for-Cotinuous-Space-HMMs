# Adaptive-Discretization-Algorithm-for-Cotinuous-Space-HMMs
his project contains functions for Poisson-AR(1) HMMs.  
My main reference is **Zimmerman et al. (2024)** [github.com/rob-zimmerman/SSM-flare](https://github.com/rob-zimmerman/SSM-flare) for some of the functions in `useful_functions.R` and `HMM_Functions.R`.  
`hmm_prediction.R` references formulas in **Zucchini et al. (2016)**.  
`auto_functions.R` is the main algorithm proposed in this project.  

## Folder Structure

- **Simulations/**  
  Contains results for the simulation studies in the report.  
  - `Experiment1/` and `Experiment3/` correspond to Experiment1 and Experiment2 in the report.

- **Earthquake/**  
  Contains codes and results for the earthquake dataset analysis, corresponding to the "Application to Earthquake Data" in the report’s Results section.  
  This folder also includes bootstrap code for parameter uncertainty evaluation.

- **MainFunctions/**  
  Contains the core functions for continuous-space HMMs, including likelihood evaluation, decoding, prediction, and the adaptive discretization algorithm.
