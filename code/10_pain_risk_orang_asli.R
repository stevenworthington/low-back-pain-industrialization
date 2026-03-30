################################################################################
#
# CAUSAL EFFECTS OF RISK FACTORS ON LOW BACK PAIN -- ORANG ASLI
# Continuous exposure effect estimation for seven risk factors
#
# This script estimates the causal effect of seven continuous risk factors on
# binary low back pain (LBP) status among the Orang Asli, using generalized
# propensity score weighting via cont_treatment_effect().
#
# Same analysis structure as 09_pain_risk_turkana.R, with key differences:
#   - Adiposity models add pregnant_ny_01 and mvpa_average_daily_min as
#     extra covariates beyond the Turkana specification
#   - CRP model adds pregnant_ny_01 and mvpa_average_daily_min
#   - PA models (steps, MVPA) add smoking, alcohol, and pregnant_ny_01
#   - Waist and body fat use energy balancing (w_meth="energy") instead
#     of optweight
#   - Adds a 7th model: occupation_subsistence_index -> LBP
#     (covariates: sex, age, bmi, smoking, alcohol, pregnant, grip_strength)
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Summary statistics for each exposure variable in its analysis sample
#   2. Model estimation -- seven continuous exposure effect models
#   3. Save results
#
# Statistical approach:
#   Covariate balancing weights (optweight or energy balancing)
#   with clarify-based simulation inference for continuous exposure effects
#   on a binary outcome.
#
# Output: models/pain_risk_orangasli.Rdata
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
# 1. SUMMARY STATISTICS FOR EACH EXPOSURE IN ITS ANALYSIS SAMPLE
# ==============================================================================

# dat_orangasli |>
#   select(low_back_pain_ny_01, crp_mg_l, sex,
#          age_years, bmi, pregnant_ny_01, smoking_ny_01,
#          alcohol_ny_01, mvpa_average_daily_min) |>
#   drop_na() |>
#   # group_by(sex) |>
#   summarize(tab = n())

# quantile(dat_orangasli[["industrialization_index"]], probs = c(0.01, 0.99), na.rm = TRUE)

# Summary of BMI restricted to its complete-case analysis sample
dat_orangasli |> 
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |> 
  select(bmi) |> 
  summary()

# Summary of waist circumference restricted to its complete-case analysis sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, waist_circum_cm, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |>
  select(waist_circum_cm) |> 
  summary()

# Summary of body fat percentage restricted to its complete-case analysis sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, body_fat_percentage, sex, age_years, pregnant_ny_01,
          mvpa_average_daily_min, market_diet_index) |>
  select(body_fat_percentage) |> 
  summary()

# Summary of log CRP restricted to its complete-case analysis sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi, pregnant_ny_01,
          smoking_ny_01, alcohol_ny_01, mvpa_average_daily_min) |>
  mutate(crp_mg_l_ln = log(crp_mg_l)) |> 
  select(crp_mg_l_ln) |> 
  summary()

# Summary of daily step counts restricted to its complete-case analysis sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, step_count_average_daily, sex, age_years, n_valid_days, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, occupation_subsistence_index, grip_strength_max_kg) |>
  select(step_count_average_daily) |> 
  summary()

# Summary of MVPA restricted to its complete-case analysis sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, occupation_subsistence_index, grip_strength_max_kg) |>
  select(mvpa_average_daily_min) |> 
  summary()

# Summary of occupation subsistence index restricted to its complete-case sample
dat_orangasli |>
  drop_na(low_back_pain_ny_01, occupation_subsistence_index, sex, age_years, bmi, smoking_ny_01,
          alcohol_ny_01, pregnant_ny_01, grip_strength_max_kg) |>
  select(occupation_subsistence_index) |>
  summary()


# ==============================================================================
# 2. MODEL ESTIMATION -- CONTINUOUS EXPOSURE EFFECT MODELS
# ==============================================================================

# BMI -> LBP: optweight, n_moments=1
# Covariates include pregnant_ny_01 and mvpa (unlike Turkana specification)
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "bmi", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "pregnant_ny_01",
                 "mvpa_average_daily_min", "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlim = c(10.70, 50.80),
  xlab = "BMI"
  ) ->
pain_risk_bmi_orangasli # n: 806, ESS: 714.1

# Waist circumference -> LBP: energy balancing (not optweight), n_moments=1
cont_treatment_effect(
  dat = dat_orangasli,
  treatment = "waist_circum_cm", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "pregnant_ny_01",
                 "mvpa_average_daily_min", "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "energy",
  xlim = c(34.00, 129.00),
  xlab = "Waist Circumference (cm)"
  ) ->
pain_risk_waist_orangasli # n: 803, ESS: 527.0

# Body fat % -> LBP: energy balancing (not optweight), n_moments=1
cont_treatment_effect(
  dat = dat_orangasli,
  treatment = "body_fat_percentage", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "pregnant_ny_01",
                 "mvpa_average_daily_min", "market_diet_index"), 
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "energy",
  xlim = c(3.4, 58.00),
  xlab = "Body Fat (%)"
  ) ->
pain_risk_fat_orangasli # n: 792, ESS: 225.0

# CRP -> LBP: log-transformed, wint=FALSE, n_moments=2
# Adds pregnant_ny_01 and mvpa_average_daily_min beyond the Turkana CRP model
cont_treatment_effect(
  dat = dat_orangasli,
  treatment = "crp_mg_l", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "bmi", "pregnant_ny_01",
                 "smoking_ny_01", "alcohol_ny_01",
                 "mvpa_average_daily_min"),
  outcome_type = "binary",
  n_moments = 2,
  wint = FALSE,
  interactions = FALSE,
  w_meth = "optweight",
  transform_expo = "log",
  xlim = c(-5.8730, 4.5770),
  xlab = "CRP (mg/L)"
  ) ->
pain_risk_crp_orangasli # n: 225, ESS: 189.6

# Daily step counts -> LBP: energy balancing, adds smoking, alcohol, pregnant beyond Turkana
cont_treatment_effect(
  dat = dat_orangasli,
  treatment = "step_count_average_daily", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "n_valid_days", "bmi",
                 "smoking_ny_01", "alcohol_ny_01",
                 "pregnant_ny_01", "occupation_subsistence_index",
                 "grip_strength_max_kg"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "energy",
  xlim = c(1084, 45753),
  xlab = "Average Daily Step Count"
  ) ->
pain_risk_steps_orangasli # n: 787, ESS: 517.7

# MVPA -> LBP: energy balancing, same expanded PA covariates
cont_treatment_effect(
  dat = dat_orangasli, 
  treatment = "mvpa_average_daily_min", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "n_valid_days", "bmi",
                 "smoking_ny_01", "alcohol_ny_01",
                 "pregnant_ny_01", "occupation_subsistence_index",
                 "grip_strength_max_kg"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "energy",
  xlim = c(0, 966.0),
  xlab = "Moderate-to-Vigorous Physical Activity"
  ) ->
pain_risk_mvpa_orangasli # n: 799, ESS: 526.3


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

# save files
save(
  pain_risk_bmi_orangasli,
  pain_risk_waist_orangasli,
  pain_risk_fat_orangasli,
  pain_risk_crp_orangasli,
  pain_risk_steps_orangasli,
  pain_risk_mvpa_orangasli,
  pain_risk_subsistence_orangasli,
  file = "models/pain_risk_orangasli.Rdata",
  compress = "gzip"
)

