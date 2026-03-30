################################################################################
#
# RISK FACTORS × INDUSTRIALIZATION — BETWEEN-GROUP DIFFERENCE
# Turkana minus Orang Asli difference in the causal effect of
# industrialization on six risk factor outcomes
#
# This script estimates the between-group difference (Turkana − Orang Asli)
# in the average dose-response function (ADRF) and average marginal effect
# function (AMEF) of industrialization on six risk factor outcomes:
#   - BMI (continuous)
#   - Waist circumference (continuous)
#   - Body fat percentage (continuous)
#   - C-reactive protein / CRP (continuous, log-transformed outcome;
#     val_seq restricted to overlapping support)
#   - Daily step counts (continuous; val_seq restricted to overlapping support)
#   - MVPA — moderate-to-vigorous physical activity (continuous;
#     val_seq restricted to overlapping support)
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Summary statistics and quantiles for industrialization index
#      (used to determine overlapping support for CRP, steps, and MVPA)
#   2. Between-group differences in causal effects for each risk factor
#   3. Save model objects
#
# Statistical approach:
#   - cont_treatment_effect_diff() fits separate optweight-balanced models
#     per group, then uses clarify simulation for difference inference
#   - n_moments=1, interactions=TRUE for all models
#   - CRP uses transform_out="log" for log-scale modeling
#   - CRP, steps, and MVPA restrict val_seq to the overlapping support
#     range of the industrialization index across both groups
#
# Dependencies: loaded in 00b_setup.R
# Helper functions: defined in 00a_functions.R
#   - cont_treatment_effect_diff(): between-group difference in ADRF/AMEF
#
# Output: models/risk_industrialization_diff.Rdata
#
################################################################################


# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================

base_dir <- "."
setwd(base_dir)

# source setup file
source("code/00b_setup.R")

set.seed(02138)


# ==============================================================================
# 1. SUMMARY STATISTICS & QUANTILES FOR INDUSTRIALIZATION INDEX
# ==============================================================================

# Summary of industrialization index for each population, restricted to
# complete cases on MVPA (the most restrictive outcome for sample overlap)

dat_turkana |>
  drop_na(mvpa_average_daily_min, industrialization_index, sex, age_years) |>
  select(industrialization_index) |>
  summary()

dat_orangasli |>
  drop_na(mvpa_average_daily_min, industrialization_index, sex, age_years) |>
  select(industrialization_index) |>
  summary()

# 1st and 99th percentiles of industrialization index among CRP complete cases
# Used to determine overlapping support range for val_seq in CRP analysis

dat_turkana |>
  drop_na(industrialization_index, crp_mg_l, sex, age_years) |>
  summarise(q01 = quantile(industrialization_index, 0.01, na.rm = TRUE),
            q99 = quantile(industrialization_index, 0.99, na.rm = TRUE))

dat_orangasli |>
  drop_na(industrialization_index, crp_mg_l, sex, age_years) |>
  summarise(q01 = quantile(industrialization_index, 0.01, na.rm = TRUE),
            q99 = quantile(industrialization_index, 0.99, na.rm = TRUE))

# ==============================================================================
# 2. BETWEEN-GROUP DIFFERENCES IN CAUSAL EFFECTS FOR EACH RISK FACTOR
# ==============================================================================

# Each block below estimates the Turkana − Orang Asli difference in ADRF/AMEF
# of industrialization on a single risk factor outcome. CRP, steps, and MVPA
# restrict the evaluation grid (val_seq) to the overlapping support range.

# BMI
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "bmi", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_bmi_diff 

# waist circumference
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "waist_circum_cm", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_waist_diff 

# percentage body fat
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "body_fat_percentage", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_fat_diff 

# CRP
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "crp_mg_l", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight",
  transform_out = "log", # log?
  val_seq = seq(8.54, 18.58, length.out=50)
  ) ->
risk_crp_diff

# daily step counts
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "step_count_average_daily", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight",
  val_seq = seq(7.365, 29.46, length.out=50)
  ) ->
risk_steps_diff

# time per day engaged in moderate-to-vigorous physical activity
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "mvpa_average_daily_min", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight",
  val_seq = seq(7.365, 29.46, length.out=50)
  ) ->
risk_mvpa_diff

# ==============================================================================
# 3. SAVE MODEL OBJECTS
# ==============================================================================

# save files
save(
  risk_bmi_diff,
  risk_waist_diff,
  risk_fat_diff,
  risk_crp_diff,
  risk_steps_diff,
  risk_mvpa_diff,
  file = "models/risk_industrialization_diff.Rdata",
  compress = "gzip"
)

