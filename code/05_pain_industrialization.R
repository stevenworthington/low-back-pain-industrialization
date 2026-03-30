
################################################################################
#
# PAIN × INDUSTRIALIZATION
# Causal effect of industrialization on low back pain (Turkana & Orang Asli)
#
# This script estimates the causal effect of industrialization (continuous
# exposure) on low back pain (binary outcome) for Turkana and Orang Asli,
# plus the between-group difference in the average dose-response function.
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Summary statistics and quantiles for industrialization index
#      (restricted to LBP complete cases in each group)
#   2. Causal effect of industrialization on LBP — Turkana
#   3. Causal effect of industrialization on LBP — Orang Asli
#   4. Between-group difference in ADRF/AMEF
#
# Statistical approach:
#   - WeightIt with optweight method for entropy balancing on sex + age
#     (n_moments=2 for mean and variance balancing)
#   - Outcome model: natural splines (df=4) for industrialization interacted
#     with covariates; Firth-type bias-reduced logistic regression (br=TRUE)
#   - G-computation via marginaleffects::avg_predictions for AERF
#   - marginaleffects::avg_slopes for AMEF
#   - Between-group ADRF/AMEF difference uses cont_treatment_effect_diff,
#     which fits separate weighted models per group and uses
#     clarify::sim + sim_adrf for difference inference with uncertainty
#
# Dependencies: loaded in 0b_setup.R
# Helper functions: defined in 0a_functions.R
#   - cont_treatment_effect(): single-group continuous exposure effect
#   - cont_treatment_effect_diff(): between-group difference in ADRF/AMEF
#
################################################################################


# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================

base_dir <- "."
setwd(base_dir)

# source setup file
source("code/0b_setup.R")

set.seed(02138)
  

# ==============================================================================
# 1. SUMMARY STATISTICS FOR INDUSTRIALIZATION INDEX
# ==============================================================================

# Summary and 1st/99th percentile quantiles for industrialization index,
# restricted to complete cases on LBP, industrialization, sex, and age

dat_turkana |>
  drop_na(low_back_pain_ny_01, industrialization_index, sex, age_years) |>
  select(industrialization_index) |>
  summary()

dat_orangasli |> 
  drop_na(low_back_pain_ny_01, industrialization_index, sex, age_years) |> 
  select(industrialization_index) |> 
  summary()

dat_turkana |>
  drop_na(low_back_pain_ny_01, industrialization_index, sex, age_years) |>
  summarise(q01 = quantile(industrialization_index, 0.01, na.rm = TRUE),
            q99 = quantile(industrialization_index, 0.99, na.rm = TRUE))                                                                       
dat_orangasli |>
  drop_na(low_back_pain_ny_01, industrialization_index, sex, age_years) |>
  summarise(q01 = quantile(industrialization_index, 0.01, na.rm = TRUE),
            q99 = quantile(industrialization_index, 0.99, na.rm = TRUE))


# ==============================================================================
# 2. CAUSAL EFFECT OF INDUSTRIALIZATION ON LBP — TURKANA
# ==============================================================================

# Estimate continuous exposure effect using optweight entropy balancing
# on sex + age, with natural spline outcome model and Firth bias reduction

# low back pain
cont_treatment_effect(
  dat = dat_turkana,
  treatment = "industrialization_index", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years"),
  outcome_type = "binary",
  n_moments = 2,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
pain_lowback_turkana # ESS: 1167, 1102.1


# save files
save(
  pain_lowback_turkana,
  file = "models/pain_industrialization_turkana.Rdata",
  compress = "gzip"
)


# ==============================================================================
# 3. CAUSAL EFFECT OF INDUSTRIALIZATION ON LBP — ORANG ASLI
# ==============================================================================

# Same approach as Section 2, applied to Orang Asli data

# low back pain
cont_treatment_effect(
  dat = dat_orangasli,
  treatment = "industrialization_index", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years"),
  outcome_type = "binary",
  n_moments = 2,
  interactions = TRUE,
  w_meth = "optweight"
) ->
pain_lowback_orangasli # ESS: 1037, 983.1


# save files
save(
  pain_lowback_orangasli,
  file = "models/pain_industrialization_orangasli.Rdata",
  compress = "gzip"
)


# ==============================================================================
# 4. BETWEEN-GROUP DIFFERENCE IN ADRF/AMEF
# ==============================================================================

# Fit separate weighted models per group and use clarify simulation
# to compute ADRF/AMEF differences with uncertainty

# low back pain
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years"),
  outcome_type = "binary",
  n_moments = 2,
  interactions = TRUE,
  w_meth = "optweight"
) ->
pain_lowback_diff 

# save files
save(
  pain_lowback_diff,
  file = "models/pain_industrialization_diff.Rdata",
  compress = "gzip"
)


