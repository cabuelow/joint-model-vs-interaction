### Simulation study to explore joint model vs. simple linear model to estimate MPA effects

This repository contains code to simulate fish biomass (in log units) and determine whether a joint modelling approach is required to estimate the effects of MPA reserve age and size, or if using an interaction (i.e., assuming 0s for MPA age and size where MPAs are not present) results is also able to accurately recover effect sizes.

- 'scripts/sim-data.R' simulates the fish biomass data, fits the different models, and plots the modelled effect size estimates
- 'script/models/joint-model.stan' is the STAN code for the joint model
