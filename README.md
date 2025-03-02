# Improved order selection method for hidden Markov models


Repository for code from the paper: "Improved order selection method for hidden Markov models: a case study with movement data". The link to the paper is available [here](https://arxiv.org/abs/2411.18826).

# Overview
The Narwhal1h2monthdts.RData file contains the tracking data used in the paper. Raw data were preprocessed to obtain one location per hour from August to late September 2017.


The rest of the repository consists of two parts:


(1) In the Simulation.R script, we provide functions to simulate data from various misspecification scenarios. We also provide details of the initialization and fitting process to select the best model using AIC and BIC criteria.

(2) The [tutorial file](Code_DPMLE_MEE/Tutorial) provides .pdf and .rmd files on how to use the DPMLE, both stationary and non-stationary, on clipped narwhal data with a single track. 



