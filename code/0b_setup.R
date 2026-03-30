################################################################################
#
# SETUP & DATA PREPARATION
# Low back pain causal analysis pipeline
#
# This script initializes the analysis environment. It is sourced at the top
# of every analysis script and performs the following steps:
#
#   1. Load all required R packages via pacman::p_load
#   2. Source shared helper functions from 0a_functions.R
#   3. Create output directories (figures/, models/) if they don't exist
#   4. Load and prepare datasets (from pre-cleaned CSVs in data/):
#      a. Turkana pain data  (clean names, convert strings to factors,
#         recode pregnancy for males, log-transform CRP, set sex factor)
#      b. Orang Asli HeLP pain data  (same pipeline as Turkana)
#      c. Combined Turkana + Orang Asli dataset  (row-bind with group column)
#      d. Tsimane back pain data  (clean names, rename columns, set sex
#         factor, add group label)
#      e. Tarahumara low back pain data  (clean names only)
#
# Dependencies:
#   pacman, tidyverse, janitor, and all packages listed in p_load() below.
#   Helper functions: 0a_functions.R
#
# Data files (in data/ directory):
#   Turkana_pain.csv
#   OA_HeLP_pain.csv
#   Tsimane_back_pain.csv
#   Tarahumara_low_back_pain.csv
#
################################################################################


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

suppressPackageStartupMessages(require(pacman))

require(pacman)
pacman::p_load(tidyverse, janitor, purrr, scales, forcats, stringr,
               patchwork, ggokabeito, cowplot, grid,
               mgcv, marginaleffects, WeightIt, cobalt, clarify)

# global options
options(scipen=20)


# ==============================================================================
# 2. UTILITY FUNCTIONS
# ==============================================================================

# Source shared function library (themes, descriptive helpers, causal inference
# functions, multi-group prevalence estimation); see 0a_functions.R for details
source("code/0a_functions.R")


# ==============================================================================
# 3. OUTPUT DIRECTORIES
# ==============================================================================

# Create output directories if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)


# ==============================================================================
# 4. DATA
# ==============================================================================

# --- 4a. Turkana pain data ---------------------------------------------------
# Load pre-cleaned Turkana dataset, convert strings to factors, 
# recode pregnancy (males = 0), log-transform CRP, set sex factor.
suppressMessages(
  suppressWarnings(
read_csv("data/Turkana_pain.csv", show_col_types=FALSE) |>
  janitor::clean_names(case="snake") |>
  mutate(across(where(is.character), factor)) |>
  mutate(
    pregnant_ny_01 = case_when(
      sex == "Male"   & is.na(pregnant_ny_01) ~ 0, # code males as zero
      sex == "Female" & pregnant_ny_01 == 1 ~ 1,
      sex == "Female" & pregnant_ny_01 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    ln_crp_mg_l = log(crp_mg_l),
    sex = factor(sex, levels = c("Female", "Male"))) ->
dat_turkana
))

# --- 4b. Orang Asli HeLP pain data ------------------------------------------
# Load pre-cleaned Orang Asli dataset, convert strings to factors, 
# recode pregnancy (males = 0), log-transform CRP, set sex factor.
suppressMessages(
  suppressWarnings(
read_csv("data/OA_HeLP_pain.csv", show_col_types=FALSE) |>
  janitor::clean_names(case="snake") |>
  mutate(across(where(is.character), factor)) |>
  mutate(
    pregnant_ny_01 = case_when(
      sex == "Male"   & is.na(pregnant_ny_01) ~ 0, # code males as zero
      sex == "Female" & pregnant_ny_01 == 1 ~ 1,
      sex == "Female" & pregnant_ny_01 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    ln_crp_mg_l = log(crp_mg_l),
    sex = factor(sex, levels = c("Female", "Male"))) ->
dat_orangasli
))

# --- 4c. Combined Turkana + Orang Asli dataset -------------------------------
# Row-bind the two datasets with a group identifier column for joint analyses.
bind_rows(
  dat_turkana   |> mutate(group = "Turkana"),
  dat_orangasli |> mutate(group = "Orang Asli")
 ) ->
dat


# --- 4d. Tsimane back pain data -----------------------------------------------
# Load pre-cleaned Tsimane dataset.
# Rename verbose column names, set sex factor levels, add group label.
read_csv("data/Tsimane_back_pain.csv", show_col_types=FALSE) |>
  janitor::clean_names(case="snake") |>
  rename(back_pain = back_pain_lasting_30_or_more_days_1_present,
         primary_cause = self_reported_primary_cause_of_back_pain) |>
  mutate(sex = factor(sex, levels = c("Female", "Male")),
         group = "Tsimane") ->
dat_tsimane

# --- 4e. Tarahumara low back pain data ----------------------------------------
# Load pre-cleaned Tarahumara dataset.
read_csv("data/Tarahumara_low_back_pain.csv", show_col_types=FALSE) |>
  janitor::clean_names(case="snake") ->
dat_tarahumara


