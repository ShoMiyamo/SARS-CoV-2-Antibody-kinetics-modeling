# ============================================================
# NT (neutralization titer) -> absolute risk lookup table
# and NT-based protection curves over time (BA.5 wave).
#
# This script:
#   1) loads a BRMS model (RDS),
#   2) builds a lookup table: anti_s (log10 NT in this model) -> absolute risk -> relative risk,
#   3) joins the lookup to time-course NT summaries by group,
#   4) plots protection (= 1 - RR) with GAM smoothing,
#   5) repeats the join/plot for multiple "NT escape" scenarios.
#
# NOTE:
# - The BRMS model RDS should be created using the Zenodo repository
#   corresponding to Miyamoto et al., Commun Med (2025).
#   Paper DOI: 10.1038/s43856-025-00894-8
#   Zenodo DOI: 10.5281/zenodo.15422823
#
# - "anti_s" in the BRMS model refers to the log10-transformed NT predictor used in that model.
#   (The name is kept for compatibility with the original modeling code.)
# ============================================================

library(dplyr)
library(brms)
library(readr)
library(ggplot2)
library(mgcv)

library(purrr)
library(tibble)

# ---------------------------
# Configuration
# ---------------------------
MODEL_RDS <- 'model/NT_GAM_symp_CommunMed.Rds' #see "Miyamoto et al., Commun Med (2025)."
NT_EST_CSV <- 'data/BA5NT_est_all_by_group.csv' #see "Hierarchical model3_infectionhistoryNT.R"

# Baseline symptomatic risk for the unvaccinated, infection-naïve group (BA.5 wave).
BASE_RISK <- 16 / (16 + 71)

# Grid resolution used in brms::conditional_effects()
CE_RESOLUTION <- 1000

# Rounding used to match predictor values between:
# - conditional_effects grid, and
# - external NT summaries (median column)
PRED_DIGITS <- 2

# Levels to keep / ordering in plots
VAX_LEVELS <- c("Twice", "3–5 times")
INF_LEVELS <- c("Uninfected", "Omicron")

# Escape scenarios: neutralization drop by a factor (log10 shift)
ESCAPE_SPECS <- tibble::tibble(
  name = c("base", "x2", "x10", "x30"),
  escape_factor = c(1, 2, 10, 30),
  title = c(
    "BA.5 NT-based protection vs symptomatic infection",
    "Predicted protection vs 2-fold NT escape variant",
    "Predicted protection vs 10-fold NT escape variant",
    "Predicted protection vs 30-fold NT escape variant"
  )
)

# Plot palette (named for stability)
PAL_VAX <- c(
  "Twice"     = "#9966b4",
  "3–5 times" = "#d3c848"
)

LT_INF <- c(
  "Uninfected" = "solid",
  "Omicron"    = "dashed"
)

# ---------------------------
# Helpers
# ---------------------------

# Find the median predictor column in the NT summary CSV.
# The original files often store the median as "X50." (legacy quantile column naming).
pick_median_col <- function(df) {
  if ("X50." %in% names(df)) return("X50.")
  if ("q50"  %in% names(df)) return("q50")
  stop("Could not find a median column. Expected 'X50.' or 'q50'.")
}

# Standardize the time column as "days".
standardize_days_col <- function(df) {
  if ("days" %in% names(df)) return(df)
  if ("X" %in% names(df))    return(df %>% rename(days = X))
  stop("Could not find a time column. Expected 'days' or legacy 'X'.")
}

# Build a lookup table from conditional_effects:
# anti_s -> absolute risk -> relative risk (RR) vs BASE_RISK.
make_risk_lookup <- function(model,
                             base_risk,
                             resolution = 1000,
                             digits = 2,
                             predictor = "anti_s") {
  
  ce <- brms::conditional_effects(
    model,
    surface    = TRUE,
    ordinary   = TRUE,
    resolution = resolution,
    ask        = FALSE,
    re_formula = NULL
  )
  
  ce_pred <- ce[[predictor]]
  if (is.null(ce_pred)) {
    stop(sprintf("conditional_effects() did not return component '%s'.", predictor))
  }
  
  # Identify the predictor value at which predicted risk is maximal.
  # We then restrict to predictor >= that value and de-duplicate by:
  # keeping, for each rounded predictor value, the row with the highest predicted risk.
  pred_at_max <- ce_pred[[predictor]][which.max(ce_pred$estimate__)]
  
  ce_pred %>%
    filter(.data[[predictor]] >= pred_at_max) %>%
    mutate(anti_s = round(.data[[predictor]], digits)) %>%  # output key is "anti_s"
    group_by(anti_s) %>%
    slice_max(estimate__, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      anti_s,
      risk_est   = estimate__,
      risk_lower = lower__,
      risk_upper = upper__,
      rr  = risk_est   / base_risk,
      rrL = risk_lower / base_risk,
      rrU = risk_upper / base_risk
    )
}

# Prepare a joined dataset for a given escape factor.
# escape_factor = 1 means no shift.
# escape_factor = 2 means a 2-fold reduction in NT -> subtract log10(2) from the median NT.
prepare_nt_protection <- function(est_df,
                                  lookup_df,
                                  escape_factor = 1,
                                  digits = 2,
                                  vax_levels = VAX_LEVELS,
                                  inf_levels = INF_LEVELS) {
  
  est_df <- standardize_days_col(est_df)
  med_col <- pick_median_col(est_df)
  
  shift <- log10(escape_factor)
  
  est_df %>%
    mutate(
      anti_s = round(.data[[med_col]] - shift, digits)
    ) %>%
    inner_join(lookup_df, by = "anti_s") %>%
    filter(infection %in% inf_levels) %>%
    mutate(
      vaccination = factor(vaccination, levels = vax_levels),
      infection   = factor(infection,   levels = inf_levels),
      protective_effect = 1 - rr
    )
}

# Plot protection over time using GAM smoothing.
plot_protection_gam <- function(df,
                                title,
                                pal_vax = PAL_VAX,
                                lt_inf = LT_INF,
                                ylim = c(0.7, 1),
                                xlim = c(0, 730)) {
  
  ggplot(
    df,
    aes(
      x = days,
      y = protective_effect,
      colour   = vaccination,
      linetype = infection,
      group    = interaction(vaccination, infection)
    )
  ) +
    geom_smooth(
      method  = "gam",
      formula = y ~ s(x, bs = "cs"),
      se      = FALSE,
      size    = 1.1,
      na.rm   = TRUE
    ) +
    scale_color_manual(values = pal_vax) +
    scale_linetype_manual(values = lt_inf) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1.0)) +
    labs(
      x = "Days post exposure",
      y = "Protective effect in BA.5 wave",
      colour   = "Vaccination",
      linetype = "Prior infection",
      title = title
    ) +
    guides(
      linetype = guide_legend(
        override.aes = list(colour = "black")
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.ticks   = element_line(colour = "black"),
      axis.text    = element_text(colour = "black", size = 12),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      plot.title   = element_text(colour = "black", size = 12)
    )
}

# ---------------------------
# 1) Load model and build lookup table
# ---------------------------
mNpntS <- readRDS(MODEL_RDS)

risk_lookup_nt <- make_risk_lookup(
  model      = mNpntS,
  base_risk  = BASE_RISK,
  resolution = CE_RESOLUTION,
  digits     = PRED_DIGITS,
  predictor  = "anti_s"
)

stopifnot(nrow(risk_lookup_nt) == dplyr::n_distinct(risk_lookup_nt$anti_s))

# ---------------------------
# 2) Load NT summaries and build datasets for each escape scenario
# ---------------------------
est_nt <- readr::read_csv(NT_EST_CSV, show_col_types = FALSE)

nt_data_list <- ESCAPE_SPECS %>%
  mutate(
    data = map(
      escape_factor,
      ~ prepare_nt_protection(
        est_df       = est_nt,
        lookup_df    = risk_lookup_nt,
        escape_factor = .x,
        digits       = PRED_DIGITS
      )
    )
  )

# Optional: keep the original-style objects if you prefer explicit names
brgNT      <- nt_data_list$data[[which(nt_data_list$name == "base")]]
brgNT_x2   <- nt_data_list$data[[which(nt_data_list$name == "x2")]]
brgNT_x10  <- nt_data_list$data[[which(nt_data_list$name == "x10")]]
brgNT_x30  <- nt_data_list$data[[which(nt_data_list$name == "x30")]]

# ---------------------------
# 3) Plot protection curves for each scenario
# ---------------------------
nt_plot_list <- nt_data_list %>%
  mutate(
    plot = map2(data, title, ~ plot_protection_gam(.x, title = .y))
  )

# Optional: keep the original-style plot objects
NTce     <- nt_plot_list$plot[[which(nt_plot_list$name == "base")]]
NTce_x2  <- nt_plot_list$plot[[which(nt_plot_list$name == "x2")]]
NTce_x10 <- nt_plot_list$plot[[which(nt_plot_list$name == "x10")]]
NTce_x30 <- nt_plot_list$plot[[which(nt_plot_list$name == "x30")]]

# Print the base plot (others are in NTce_x2 / NTce_x10 / NTce_x30)
NTce
NTce_x2
NTce_x10
NTce_x30