## R/Stan codes for SARS-CoV-2 Antibody kinetics modeling
[![DOI](https://zenodo.org/badge/1147695304.svg)](https://doi.org/10.5281/zenodo.18454709)

This repository provides R/Stan scripts to model longitudinal SARS‑CoV‑2 antibody kinetics
(e.g., anti‑Spike binding antibody and neutralizing titers) and to translate antibody levels
into predicted protection against symptomatic infection.

Please cite the associated preprint when using this code:
- [Numakura & Miyamoto et al., medRxiv, 2025](https://doi.org/10.1101/2025.11.14.25340279)

Some downstream “correlates of protection” steps use pre‑fitted brms/GAM models
from Miyamoto et al., Communications Medicine (2025):
- [Miyamoto et al., Commun Med, 2025](https://doi.org/10.1038/s43856-025-00894-8)


## What this repository does

1. **Antibody kinetics (Bayesian models in Stan)**
   - Fits exponential‑decay–type models to log10-transformed antibody titers.
   - Generates posterior summaries and predictive draws.
   - Supports stratification by immune history (vaccination / prior infection), age/sex groups,
     and bivalent booster settings (depending on the script).

2. **Risk / protection mapping (anti‑S or NT → absolute risk / relative risk)**
   - Uses `brms::conditional_effects()` outputs from pre‑fitted GAM models to build a
     lookup table: antibody level ↔ predicted symptomatic risk.
   - Joins this lookup with the time‑trajectory estimates from the Stan models to plot
     “protection vs time”.
   - Includes optional “immune escape” scenarios for NT by applying log10 fold‑reductions
     (e.g., 2‑fold, 10‑fold, 30‑fold reductions).


## Repository structure

- `stan/`  
  Stan model files (e.g., exponential decay and group‑structured variants).

- `script1/`  
  Antibody kinetics workflow (Stan fitting + post‑processing + figures/tables).
  Scripts here typically:
  - read the cleaned input data,
  - fit a Stan model (`rstan::sampling()`),
  - extract posterior draws (`rstan::extract()`),
  - create quantile summary tables used for plotting and export.

- `script2/`  
  Protection / risk workflow (anti‑S or NT → predicted risk/protection + plots).
  Scripts here typically:
  - load a pre‑fitted `brms` model object (`.Rds`),
  - obtain conditional effects surfaces,
  - build antibody ↔ risk lookup tables,
  - merge with kinetics summaries and plot protection vs time.

- `data/`  
  Input data used by the scripts (CSV).  
  
- `model/`  
  Optional: saved fitted objects (`.Rds`) or intermediate outputs used by scripts.

## Data requirements 

**Source data (antibody kinetics)**
The source dataset used for antibody-kinetics fitting is provided as Supplementary/Source Data in the manuscript.


## Dependencies
- R version 4.5.2 (2025-10-31)

- R packages:
  - rstan_2.32.7
  - readr_2.1.5
  - ggplot2_4.0.1
  - dplyr_1.1.4
  - cowplot_1.2.0
  - purrr_1.1.0
  - tibble_3.3.0
  - brms_2.23.0
  - mgcv_1.9-3
