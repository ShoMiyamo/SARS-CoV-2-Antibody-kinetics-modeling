# R codes for SARS-CoV-2-Antibody-kinetics-modeling

## Summary
This repository contains computational codes used in [Miyamoto et al., medRxiv, 2025](https://doi.org/10.1101/2025.11.14.25340279).

## Contents
- script
  - Anti_N antibody_response.R : Modeling the antibody response main script
    - nmodel01.stan : Modeling the antibody response stan model file
  - infection_risk_estimation.R


## Dependencies
- R version 4.5.2 (2025-10-31)

- R packages:
  - rstan_2.32.7
  - readr_2.1.5
  - ggplot2_4.0.1
  - dplyr_1.1.4
  - cowplot_1.2.0
