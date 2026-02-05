# Navigating Illusions of Importance

This repository contains code and data sufficient to reproduce the simulation experiments and the empirical river flow analysis presented in the manuscript
"Interpreting Random Forests and Navigating Illusions of Importance in Spatial Data."

The repository is organized into two main components:

- SimStudy/: controlled simulation experiments
- RiverSeasonality/: empirical analysis of river flow seasonality across North American catchments

---

## Repository structure

### 1. SimStudy/ - Simulation experiments

The SimStudy/ directory contains code for the linear and nonlinear simulation studies used to examine how spatial dependence affects random forest interpretation.

Main scripts:
- Linear_Sim.R  
  Runs the linear simulation study with spatially structured and unstructured covariates.

- Nonlinear_Sim.R  
  Runs the nonlinear and mixed-type simulation study.

Helper functions:
- fitMaternGP.R  
  Functions for fitting Gaussian process models with a Matern covariance, used to estimate spatial correlation parameters.

- mkNNIndx.R  
  Utility function for constructing nearest-neighbor index sets for Vecchia-based spatial approximations.

- TransformFunctions.R  
  Functions implementing spatial whitening (decorrelation) and the corresponding back-transformation for prediction.

Each simulation script generates a single representative simulation replicate, sufficient to demonstrate differences in predictive performance and variable importance between models that ignore spatial dependence and models that account for it.

---

### 2. RiverSeasonality/ - Empirical river flow analysis

The RiverSeasonality/ directory contains the data and code used for the empirical analysis of river flow seasonality across North American catchments.

Files:
- PCA.RData  
  Serialized R data file containing the river seasonality index and associated catchment-level covariates used in the analysis.

- river_vip_comparison_parallel.R  
  Script that fits independent and spatially decorrelated random forest models, evaluates predictive performance, and computes variable importance comparisons reported in the manuscript.

This analysis reproduces the river flow application results, including model performance metrics and variable importance patterns comparing spatial and non-spatial random forests.

---

## How to run the analyses

Each component can be run independently.

Simulation studies:
1. Set the working directory to SimStudy/
2. Run either:
   - Linear_Sim.R, or
   - Nonlinear_Sim.R

River seasonality analysis:
1. Set the working directory to RiverSeasonality/
2. Load the data and run:
   - river_vip_comparison_parallel.R

---

## Software requirements

All analyses were developed and run using:

- R (version 4.4.0 or higher)

Key R packages include:
- ranger
- Matrix
- fields
- stats
- parallel

Because some analyses involve random number generation, numerical results may vary slightly unless a fixed random seed is used (seeds are set within scripts where appropriate).

---

## Data

- All data used in the simulation studies are fully simulated.
- The river seasonality analysis uses derived catchment-level data stored in PCA.RData, which were constructed from publicly available hydrologic and environmental data sources cited in the manuscript.

No restricted or proprietary datasets are included in this repository.

---


## License

This repository is released under an open-source license for reproducibility and reuse.