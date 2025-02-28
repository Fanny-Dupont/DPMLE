# Tutorial

This repository contains a tutorial for analyzing narwhal movement data using both DPMLE methods. The tutorial focuses on fitting stationary and non-stationary models, incorporating covariates, and interpreting model outputs.

## Files Description

- **OneNarwhal_Tutorial.RData**:
  - This file contains the dataset used for the tutorial. It consists of a clipped track of narwhal movement data with the identifier track `ID=10`.
  
- **Tutorial.Rmd**:
  - This R Markdown file provides the code and narrative for the tutorial. It provides a step-by-step guide on fitting models, incorporating covariates, and interpreting results using the narwhal data.

- **Tutorial.pdf**:
  - This is the compiled PDF version of the `Tutorial.Rmd` file. It presents the tutorial in a readable format, suitable for distribution or reference.

- **sourcefunctions.R**:
  - This file contains the source functions required for the tutorial. These functions are used in the `Tutorial.Rmd` file to perform various analyses and model fitting tasks. It contains all the functions required to fit the EM algorithm for gamma and von Mises state-dependent distributions (e.g., forward algorithm, (I) and (II) in eq (14) of the paper).

- **tutorial_source.R**:
  - This file includes additional source codes that support the tutorial. It contains the EM algorithm to fit both DPMLE methods with gamma and von mises state-dependent distributions. It can be easily adapted to any other state-dependent distributions.

## Tutorial Learning Objectives

1. Fit stationary DPMLE on movement data.
2. Interpret the model output and extract estimates.
3. Check model fit.
4. Incorporate covariates on behaviour transition probabilities with non-stationary DPMLE.
5. Interpret the model output and extract estimates.
6. Check model fit.
7. Compare both models, select the best performing one.
8. (Bonus) Select the best tuning parameters (explore 10) with BIC_DPMLE for stationary DPMLE.

## Getting Started

To get started with the tutorial, follow these steps:

1. Clone the repository to your local machine.
2. Open the `Tutorial.Rmd` file in RStudio.
3. Ensure you have the necessary R packages installed.
4. Knit the R Markdown file to generate the tutorial PDF.
5. Follow the instructions in the tutorial to analyze the narwhal movement data.

## Requirements

- R and RStudio
- Necessary R packages (listed in the `Tutorial.Rmd` file)

## Contact
For any questions or issues, please contact Fanny Dupont at fanny.dupont@stat.ubc.ca.

