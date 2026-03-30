################################################################################
#
# CAUSAL EFFECTS OF RISK FACTORS ON LOW BACK PAIN -- TURKANA
# Continuous exposure effect estimation for six risk factors
#
# This script estimates the causal effect of six continuous risk factors on
# binary low back pain (LBP) status among the Turkana, using generalized
# propensity score weighting via cont_treatment_effect().
#
# Each risk factor uses a tailored set of confounders reflecting the assumed
# causal structure for that specific exposure-outcome pair:
#   - Adiposity measures (BMI, waist, body fat): sex, age, market_diet_index
#   - CRP: sex, age, bmi, smoking, alcohol (log-transformed exposure)
#   - Physical activity (steps, MVPA): sex, age, n_valid_days, bmi, alcohol,
#     occupation_subsistence_index, grip_strength
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Data exploration -- complete case counts, missing data, summary stats
#   2. Summary statistics for each exposure variable in its analysis sample
#   3. Model estimation -- six continuous exposure effect models
#   4. Save results
#
# Statistical approach:
#   Covariate balancing weights (optweight or energy balancing)
#   with clarify-based simulation inference for continuous exposure effects
#   on a binary outcome. Balancing weights target covariate moments (n_moments)
#   to achieve approximate independence of exposure and confounders.
#
# Output: models/pain_risk_turkana.Rdata
# Dependencies: loaded in 0b_setup.R
# Helper functions: defined in 0a_functions.R (cont_treatment_effect)
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
# 1. DATA EXPLORATION -- COMPLETE CASE COUNTS & MISSING DATA
# ==============================================================================

# Complete case count for all risk factor variables combined
dat_turkana |>
  select(low_back_pain_ny_01, crp_mg_l, sex,
         age_years, bmi, pregnant_ny_01, smoking_ny_01,
         alcohol_ny_01, mvpa_average_daily_min,
         occupation_subsistence_index, market_diet_index,
         grip_strength_max_kg, n_valid_days) |>
  drop_na() |>
  # group_by(sex) |>
  summarize(tab = n())

# quantile(dat_turkana[["industrialization_index"]], probs = c(0.01, 0.99), na.rm = TRUE)

# Non-missing counts per variable
dat_turkana |>
  select(low_back_pain_ny_01, crp_mg_l, sex,
         age_years, bmi, pregnant_ny_01, smoking_ny_01,
         alcohol_ny_01, mvpa_average_daily_min,
         occupation_subsistence_index, market_diet_index,
         grip_strength_max_kg, n_valid_days) |>
  summarize(across(everything(), ~ sum(!is.na(.)))) |> 
  pivot_longer(cols = everything())


# Complete cases for adiposity model covariates
dat_turkana |>
  select(low_back_pain_ny_01, bmi, sex,
         age_years, pregnant_ny_01, market_diet_index) |>
  drop_na() |>
  #group_by(sex) |>
  summarize(tab = n())

# Complete cases for physical activity model covariates
dat_turkana |>
  select(low_back_pain_ny_01, bmi, sex,
         age_years, n_valid_days,
         smoking_ny_01, alcohol_ny_01,
         occupation_subsistence_index,
         grip_strength_max_kg
         ) |>
  drop_na() |>
  #group_by(sex) |>
  summarize(tab = n()) # n=201

# Count of smokers among participants with valid accelerometer data
dat_turkana |>
  drop_na(n_valid_days) |>
  select(smoking_ny_01) |>
  drop_na() |>
  count(smoking_ny_01)


# ==============================================================================
# 2. SUMMARY STATISTICS FOR EACH EXPOSURE IN ITS ANALYSIS SAMPLE
# ==============================================================================

# Summary of BMI restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |> 
  select(bmi) |> 
  summary()

# Summary of waist circumference restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, waist_circum_cm, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |> 
  select(waist_circum_cm) |> 
  summary()

# Summary of body fat percentage restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, body_fat_percentage, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |> 
  select(body_fat_percentage) |> 
  summary()

# Summary of log CRP restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi, pregnant_ny_01,
          smoking_ny_01, alcohol_ny_01, mvpa_average_daily_min) |>
  mutate(crp_mg_l_ln = log(crp_mg_l)) |> 
  select(crp_mg_l_ln) |> 
  summary()

# Summary of daily step counts restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, step_count_average_daily, sex, age_years, n_valid_days, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, occupation_subsistence_index, grip_strength_max_kg) |>
  select(step_count_average_daily) |>
  summary()

# Summary of MVPA restricted to its complete-case analysis sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, occupation_subsistence_index, grip_strength_max_kg) |>
  select(mvpa_average_daily_min) |>
  summary()

# Summary of occupation subsistence index restricted to its complete-case sample
dat_turkana |>
  drop_na(low_back_pain_ny_01, occupation_subsistence_index, sex, age_years, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, grip_strength_max_kg) |>
  select(occupation_subsistence_index) |>
  summary()

# ==============================================================================
# 3. MODEL ESTIMATION -- CONTINUOUS EXPOSURE EFFECT MODELS
# ==============================================================================

# BMI -> LBP: optweight, n_moments=1, no interactions
cont_treatment_effect(
  dat = dat_turkana, 
  treatment = "bmi", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", 
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlim = c(10.70, 50.80), 
  xlab = "BMI"
  ) ->
pain_risk_bmi_turkana # n: 1127, ESS: 1015.4

# Waist circumference -> LBP: optweight, n_moments=2
cont_treatment_effect(
  dat = dat_turkana, 
  treatment = "waist_circum_cm", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", 
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 2,
  interactions = FALSE,
  w_meth = "optweight",
  xlim = c(34.00, 129.00),
  xlab = "Waist Circumference (cm)"
  ) ->
pain_risk_waist_turkana # n: 1119, ESS: 1074.8

# Body fat % -> LBP: optweight, n_moments=2
cont_treatment_effect(
  dat = dat_turkana, 
  treatment = "body_fat_percentage", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", 
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 2,
  interactions = FALSE,
  w_meth = "optweight",
  xlim = c(3.4, 58.00),
  xlab = "Body Fat (%)"
  ) ->
pain_risk_fat_turkana # n: 1010, ESS: 668.2

# CRP -> LBP: log-transformed exposure, wint=FALSE, n_moments=2
# Different covariate set reflecting CRP-specific confounders
cont_treatment_effect(
  dat = dat_turkana,
  treatment = "crp_mg_l", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "bmi", 
                 "smoking_ny_01", "alcohol_ny_01"),
  outcome_type = "binary",
  n_moments = 2,
  wint = FALSE,
  interactions = FALSE,
  w_meth = "optweight",
  transform_expo = "log",
  xlim = c(-5.8730, 4.5770),
  xlab = "CRP (mg/L)"
  ) ->
pain_risk_crp_turkana # n: 361, ESS: 336.8

# Daily step counts -> LBP: PA-specific covariates including accelerometer wear days
cont_treatment_effect(
  dat = dat_turkana, 
  treatment = "step_count_average_daily", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "n_valid_days", "bmi",
                 "alcohol_ny_01",
                 "occupation_subsistence_index",
                 "grip_strength_max_kg"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlim = c(1084, 45753),
  xlab = "Average Daily Step Count"
  ) ->
pain_risk_steps_turkana # n: 201, ESS: 122.0

# MVPA -> LBP: energy balancing method (w_meth="energy"), same PA covariates
cont_treatment_effect(
  dat = dat_turkana, 
  treatment = "mvpa_average_daily_min", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "n_valid_days", "bmi",
                 "alcohol_ny_01",
                 "occupation_subsistence_index",
                 "grip_strength_max_kg"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "energy",
  xlim = c(0, 966.0),
  xlab = "Moderate-to-Vigorous Physical Activity"
  ) ->
pain_risk_mvpa_turkana # n: 201, ESS: 109.3

# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

# save files
save(
  pain_risk_bmi_turkana,
  pain_risk_waist_turkana,
  pain_risk_fat_turkana,
  pain_risk_crp_turkana,
  pain_risk_steps_turkana,
  pain_risk_mvpa_turkana,
  file = "models/pain_risk_turkana.Rdata",
  compress = "gzip"
)

