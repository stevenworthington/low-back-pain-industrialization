################################################################################
#
# RISK FACTORS × INDUSTRIALIZATION — ORANG ASLI
# Causal effect of industrialization on six risk factor outcomes (Orang Asli)
#
# This script estimates the causal effect of industrialization (continuous
# exposure) on six risk factor outcomes among the Orang Asli:
#   - BMI (continuous)
#   - Waist circumference (continuous)
#   - Body fat percentage (continuous)
#   - C-reactive protein / CRP (continuous, log-transformed outcome)
#   - Daily step counts (continuous, with n_valid_days as extra covariate)
#   - MVPA — moderate-to-vigorous physical activity (continuous, with
#     n_valid_days as extra covariate)
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Summary statistics for industrialization index
#   2. Causal effects of industrialization on each risk factor
#   3. Save model objects
#
# Statistical approach:
#   - cont_treatment_effect() with optweight entropy balancing on
#     sex + age_years (steps and MVPA additionally balance n_valid_days)
#   - n_moments=1, interactions=TRUE for all models
#   - CRP uses transform_out="log" for log-scale modeling
#   - Effective sample sizes (ESS) are noted inline for each model
#
# Dependencies: loaded in 00b_setup.R
# Helper functions: defined in 00a_functions.R
#   - cont_treatment_effect(): single-group continuous exposure effect
#
# Output: models/risk_industrialization_orangasli.Rdata
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
# 1. SUMMARY STATISTICS FOR INDUSTRIALIZATION INDEX
# ==============================================================================

# Summary of industrialization index for each population (restricted to
# complete cases on industrialization, sex, and age)

dat_orangasli |>
  drop_na(industrialization_index, sex, age_years) |>
  select(industrialization_index) |>
  summary()

dat_turkana |>
  drop_na(industrialization_index, sex, age_years) |>
  select(industrialization_index) |>
  summary()


# ==============================================================================
# 2. CAUSAL EFFECTS OF INDUSTRIALIZATION ON EACH RISK FACTOR
# ==============================================================================

# Each block below estimates the average dose-response function (ADRF) and
# average marginal effect function (AMEF) of industrialization on a single
# risk factor outcome, using optweight balancing weights.

# BMI
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "bmi", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight",
  ) ->
risk_bmi_orangasli # n: 1037, ESS: 991.3

# waist circumference
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "waist_circum_cm", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_waist_orangasli # n: 1032, ESS: 983.0

# percentage body fat
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "body_fat_percentage", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_fat_orangasli # n: 1020, ESS: 972.8

# CRP
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "crp_mg_l", 
  covariates = c("sex", "age_years"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight",
  transform_out = "log" 
  ) ->
risk_crp_orangasli # n: 237, ESS: 206.6

# daily step counts
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "step_count_average_daily", 
  covariates = c("sex", "age_years", "n_valid_days"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_steps_orangasli # n: 860, ESS: 790.5

# time per day engaged in moderate-to-vigorous physical activity
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "industrialization_index", 
  outcome = "mvpa_average_daily_min", 
  covariates = c("sex", "age_years", "n_valid_days"),
  outcome_type = "continuous",
  n_moments = 1,
  interactions = TRUE,
  w_meth = "optweight"
  ) ->
risk_mvpa_orangasli # n: 872, ESS: 803.3

# ==============================================================================
# 3. SAVE MODEL OBJECTS
# ==============================================================================

# save files
save(
  risk_bmi_orangasli,
  risk_waist_orangasli,
  risk_fat_orangasli,
  risk_crp_orangasli,
  risk_steps_orangasli,
  risk_mvpa_orangasli,
  file = "models/risk_industrialization_orangasli.Rdata",
  compress = "gzip"
)

