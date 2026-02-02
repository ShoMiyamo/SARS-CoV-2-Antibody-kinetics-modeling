# ============================================================
# Build an "anti-S ↔ absolute risk" lookup table
# and plot anti-S-based protection over time (BA.5 wave).
#
# NOTE:
# - The BRMS model RDS should be created using the Zenodo repository
#   corresponding to Miyamoto et al., Commun Med (2025).
#   Paper DOI: 10.1038/s43856-025-00894-8
#   Zenodo DOI: 10.5281/zenodo.15422823
# ============================================================

library(dplyr)
library(readr)
library(brms)
library(ggplot2)

# ---------------------------
# User configuration
# ---------------------------
BRMS_MODEL_RDS <- "model/antiS_GAM_symp_CommunMed.Rds" #see "Miyamoto et al., Commun Med (2025)."
ANTI_S_EST_CSV <- "data/anti_s_est_all_by_group.csv" #see "Hierarchical model3_infectionhistory.R"

# Baseline symptomatic risk in the BA.5 wave for the unvaccinated, infection-naïve group.
BASE_RISK <- 16 / (16 + 71)

# Rounding used to match anti_s between the lookup grid and the external summary table.
ANTI_S_DIGITS <- 2

# Grid resolution used in brms::conditional_effects().
CE_RESOLUTION <- 1000

# ---------------------------
# Helper: anti_s -> risk lookup
# ---------------------------
make_risk_lookup <- function(model,
                             base_risk,
                             resolution = 1000,
                             digits = 2) {
  # Compute conditional effects for anti_s.
  # If the model includes additional predictors, the returned object may contain a grid.
  # We construct a one-to-one mapping by:
  #   1) restricting to the region beyond the anti_s value at which predicted risk is maximal, and
  #   2) for each rounded anti_s, keeping the row with the highest predicted risk
  #      (a conservative "worst-case" mapping across the grid).
  ce <- brms::conditional_effects(
    model,
    surface    = TRUE,
    ordinary   = TRUE,
    resolution = resolution,
    ask        = FALSE,
    re_formula = NULL
  )
  
  ce_s <- ce[["anti_s"]]
  if (is.null(ce_s)) stop("conditional_effects() did not return an 'anti_s' component.")
  
  anti_s_at_max_risk <- ce_s$anti_s[which.max(ce_s$estimate__)]
  
  ce_s %>%
    filter(anti_s >= anti_s_at_max_risk) %>%
    mutate(anti_s = round(anti_s, digits)) %>%
    group_by(anti_s) %>%
    slice_max(estimate__, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      anti_s,
      risk_est   = estimate__,
      risk_lower = lower__,
      risk_upper = upper__,
      rr         = risk_est / base_risk
    )
}

# ---------------------------
# 1) Load BRMS model and build lookup table
# ---------------------------
mNs <- readRDS(BRMS_MODEL_RDS)

risk_lookup <- make_risk_lookup(
  model      = mNs,
  base_risk  = BASE_RISK,
  resolution = CE_RESOLUTION,
  digits     = ANTI_S_DIGITS
)

# Sanity check: anti_s should be unique after de-duplication
stopifnot(nrow(risk_lookup) == n_distinct(risk_lookup$anti_s))

# ---------------------------
# 2) Load anti_s time-course summaries and join with lookup
# ---------------------------
est <- readr::read_csv(ANTI_S_EST_CSV, show_col_types = FALSE)

# Accept both "days" or legacy "X" as the time variable
if ("X" %in% names(est) && !("days" %in% names(est))) {
  est <- est %>% rename(days = X)
}

# Match anti_s values by rounding the median estimate (q50)
est <- est %>%
  mutate(anti_s = round(q50, ANTI_S_DIGITS))

# Join: keep only rows where anti_s exists in the risk lookup
brg <- est %>%
  inner_join(risk_lookup, by = "anti_s") %>%
  filter(infection %in% c("Uninfected", "Omicron")) %>%
  mutate(
    vaccination = factor(vaccination, levels = c("Twice", "3 times", "4–5 times")),
    infection   = factor(infection,   levels = c("Uninfected", "Omicron")),
    protective_effect = 1 - rr
  )

# ---------------------------
# 3) Plot: protection vs days (GAM smoothing)
# ---------------------------
pal_vax <- c(
  "Twice"     = "#9966b4",
  "3 times"   = "#39568C",
  "4–5 times" = "#d3c848"
)

linetype_inf <- c(
  "Uninfected" = "solid",
  "Omicron"    = "dashed"
)

p <- ggplot(
  brg,
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
  scale_linetype_manual(values = linetype_inf) +
  coord_cartesian(ylim = c(0.7, 1)) +
  scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1.0)) +
  labs(
    x = "Days post exposure",
    y = "Protective effect in BA.5 wave",
    colour   = "Vaccination",
    linetype = "Prior infection",
    title = "Anti-S-based protection vs symptomatic infection"
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(colour = "black")
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.ticks  = element_line(colour = "black"),
    axis.text   = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title  = element_text(colour = "black", size = 12)
  )

p
