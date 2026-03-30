################################################################################
#
# BETWEEN-GROUP DIFFERENCES IN CAUSAL EFFECTS OF RISK FACTORS ON LBP
# Turkana minus Orang Asli difference estimation for six risk factors
#
# This script estimates the between-group difference (Turkana - Orang Asli)
# in the causal effect of six continuous risk factors on binary low back pain,
# using cont_treatment_effect_diff() with clarify simulation.
#
# NOTE: only one female was pregnant out of a sample of 1127, which is why
# pregnancy is not included as a covariate in the shared models.
#
# Covariates match the Turkana specification (simpler set) since they must
# be shared across both groups for valid comparison.
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Data exploration -- summary statistics for each exposure in both groups
#   2. Quantile computation -- 1st/99th percentiles for each exposure in
#      both groups to determine overlapping val_seq ranges
#   3. Model estimation -- six between-group difference models
#   4. Save results
#
# Key convergence notes:
#   - CRP: log-transformed exposure, wint=FALSE required for convergence
#   - Step counts & MVPA: smoking_ny_01 excluded from covariates due to
#     non-convergence (only 11 smokers in the combined sample)
#
# Statistical approach:
#   For each risk factor, separate balancing-weight models are fit within
#   each group, then clarify simulation computes the difference
#   in dose-response curves over a shared val_seq range (overlap of both
#   groups' 1st-99th percentile ranges).
#
# Output: models/pain_risk_diff.Rdata
# Dependencies: loaded in 00b_setup.R
# Helper functions: defined in 00a_functions.R (cont_treatment_effect_diff)
#
################################################################################


# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================

# NOTE: only one female was pregnant out of a sample of 1127

base_dir <- "."
setwd(base_dir)

# source setup file
source("00b_setup.R")

set.seed(02138)


# ==============================================================================
# 1. DATA EXPLORATION -- SUMMARY STATISTICS FOR EACH EXPOSURE IN BOTH GROUPS
# ==============================================================================

# BMI summaries for both groups (restricted to shared covariate complete cases)
dat_turkana |>
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, market_diet_index) |>
  select(bmi) |>
  summary()

dat_orangasli |>
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, market_diet_index) |>
  select(bmi) |>
  summary()

# CRP summaries (log-transformed) for both groups
dat_turkana |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi, smoking_ny_01, alcohol_ny_01) |> 
  select(crp_mg_l) |> 
  mutate(crp_mg_l = log(crp_mg_l)) |>
  summary()

dat_orangasli |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi, smoking_ny_01, alcohol_ny_01) |>
  select(crp_mg_l) |>
  mutate(crp_mg_l = log(crp_mg_l)) |>
  summary()

# MVPA summaries for both groups
dat_turkana |>
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi, alcohol_ny_01, occupation_subsistence_index, grip_strength_max_kg) |> 
  select(mvpa_average_daily_min) |> 
  summary()

dat_orangasli |> 
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi, alcohol_ny_01, occupation_subsistence_index, grip_strength_max_kg) |> 
  select(mvpa_average_daily_min) |> 
  summary()


# ==============================================================================
# 2. QUANTILE COMPUTATION -- DETERMINE OVERLAPPING VAL_SEQ RANGES
# ==============================================================================
# Compute 1st and 99th percentiles of each exposure in both groups.
# The overlapping range defines the val_seq for between-group comparison,
# ensuring both groups have adequate data support across the evaluation range.

# BMI quantiles
dat_turkana |>
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(bmi, 0.01, na.rm = TRUE),
            q99 = quantile(bmi, 0.99, na.rm = TRUE))                                                                       
dat_orangasli |>
  drop_na(low_back_pain_ny_01, bmi, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(bmi, 0.01, na.rm = TRUE),
            q99 = quantile(bmi, 0.99, na.rm = TRUE))

# Waist circumference quantiles
dat_turkana |>
  drop_na(low_back_pain_ny_01, waist_circum_cm, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(waist_circum_cm, 0.01, na.rm = TRUE),
            q99 = quantile(waist_circum_cm, 0.99, na.rm = TRUE))                                                                       
dat_orangasli |>
  drop_na(low_back_pain_ny_01, waist_circum_cm, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(waist_circum_cm, 0.01, na.rm = TRUE),
            q99 = quantile(waist_circum_cm, 0.99, na.rm = TRUE))

# Body fat percentage quantiles
dat_turkana |>
  drop_na(low_back_pain_ny_01, body_fat_percentage, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(body_fat_percentage, 0.01, na.rm = TRUE),
            q99 = quantile(body_fat_percentage, 0.99, na.rm = TRUE))                                                                       
dat_orangasli |>
  drop_na(low_back_pain_ny_01, body_fat_percentage, sex, age_years, market_diet_index) |>
  summarise(q01 = quantile(body_fat_percentage, 0.01, na.rm = TRUE),
            q99 = quantile(body_fat_percentage, 0.99, na.rm = TRUE))

# CRP quantiles (log scale)
dat_turkana |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi,
          smoking_ny_01, alcohol_ny_01) |>
  summarise(q01 = quantile(log(crp_mg_l), 0.01, na.rm = TRUE),
            q99 = quantile(log(crp_mg_l), 0.99, na.rm = TRUE)) 

dat_orangasli |>
  drop_na(low_back_pain_ny_01, crp_mg_l, sex, age_years, bmi,
          smoking_ny_01, alcohol_ny_01) |>
  summarise(q01 = quantile(log(crp_mg_l), 0.01, na.rm = TRUE),
            q99 = quantile(log(crp_mg_l), 0.99, na.rm = TRUE))

# Daily step count quantiles
dat_turkana |>
  drop_na(low_back_pain_ny_01, step_count_average_daily, sex, age_years, n_valid_days, bmi,
          alcohol_ny_01, occupation_subsistence_index,
          grip_strength_max_kg) |>
  summarise(q01 = quantile(step_count_average_daily, 0.01, na.rm = TRUE),
            q99 = quantile(step_count_average_daily, 0.99, na.rm = TRUE)) 

dat_orangasli |>
  drop_na(low_back_pain_ny_01, step_count_average_daily, sex, age_years, n_valid_days, bmi,
          alcohol_ny_01, occupation_subsistence_index,
          grip_strength_max_kg) |>
  summarise(q01 = quantile(step_count_average_daily, 0.01, na.rm = TRUE),
            q99 = quantile(step_count_average_daily, 0.99, na.rm = TRUE))

# MVPA quantiles
dat_turkana |>
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi,
          alcohol_ny_01, occupation_subsistence_index,
          grip_strength_max_kg) |>
  summarise(q01 = quantile(mvpa_average_daily_min, 0.01, na.rm = TRUE),
            q99 = quantile(mvpa_average_daily_min, 0.99, na.rm = TRUE)) 

dat_orangasli |>
  drop_na(low_back_pain_ny_01, mvpa_average_daily_min, sex, age_years, n_valid_days, bmi,
          alcohol_ny_01, occupation_subsistence_index,
          grip_strength_max_kg) |>
  summarise(q01 = quantile(mvpa_average_daily_min, 0.01, na.rm = TRUE),
            q99 = quantile(mvpa_average_daily_min, 0.99, na.rm = TRUE))


# ==============================================================================
# 3. MODEL ESTIMATION -- BETWEEN-GROUP DIFFERENCE MODELS
# ==============================================================================
# Each model fits separate balancing-weight models within Turkana (dat1) and Orang Asli
# (dat2), then uses clarify simulation to estimate the difference in
# dose-response curves. Covariates use the Turkana (simpler) specification
# to enable valid between-group comparison.

# BMI: Turkana - Orang Asli difference in BMI -> LBP effect
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "bmi", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years",
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlab = "BMI",
  xlim = c(10.70, 50.80),
  val_seq = seq(14.30, 31.30, length.out=50)
  ) ->
pain_risk_bmi_diff 

# Waist circumference: Turkana - Orang Asli difference
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "waist_circum_cm", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years",
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlab = "Waist Circumference",
  xlim = c(34.00, 132.00),
  val_seq = seq(58.4, 105, length.out=50)
  ) ->
pain_risk_waist_diff 

# Body fat %: Turkana - Orang Asli difference
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "body_fat_percentage", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years",
                 "market_diet_index"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight",
  xlab = "Body Fat %",
  xlim = c(3.4, 58.00),
  val_seq = seq(5.31, 48.3, length.out=50)
  ) ->
pain_risk_fat_diff 

# CRP: log-transformed, wint=FALSE required for convergence
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "crp_mg_l", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "bmi", 
                 "smoking_ny_01", "alcohol_ny_01"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  wint = FALSE, 
  w_meth = "optweight",
  transform_expo = "log",
  xlab = "CRP (mg/L)",
  xlim = c(-6.5966, 5.1096),
  val_seq = seq(-0.693, 3.43, length.out=50)
  ) ->
pain_risk_crp_diff  

# Daily step counts: smoking_ny_01 excluded -- causes non-convergence (only 11 smokers)
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
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
  xlab = "Average Daily Step Count",
  xlim = c(1084, 45753),
  val_seq = seq(2270, 22231, length.out=50)
  ) ->
pain_risk_steps_diff

# MVPA: smoking_ny_01 excluded -- causes non-convergence (only 11 smokers)
cont_treatment_effect_diff(
  dat1 = dat_turkana,
  dat2 = dat_orangasli, 
  treatment = "mvpa_average_daily_min", 
  outcome = "low_back_pain_ny_01", 
  covariates = c("sex", "age_years", "n_valid_days", "bmi",
                 "alcohol_ny_01",
                 "occupation_subsistence_index",
                 "grip_strength_max_kg"),
  outcome_type = "binary",
  n_moments = 1,
  interactions = FALSE,
  w_meth = "optweight", 
  xlab = "MVPA Daily Minimum",
  xlim = c(0, 966.0),
  val_seq = seq(17.4, 276, length.out=50)
  ) ->
pain_risk_mvpa_diff

# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

# save files
save(
  pain_risk_bmi_diff,
  pain_risk_waist_diff,
  pain_risk_fat_diff,
  pain_risk_crp_diff,
  pain_risk_steps_diff,
  pain_risk_mvpa_diff,
  file = "models/pain_risk_diff.Rdata",
  compress = "gzip"
)

