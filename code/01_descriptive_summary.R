
################################################################################
#
# DESCRIPTIVE SUMMARY
# Basic descriptive statistics for Turkana, Orang Asli, Tsimane, & Tarahumara
#
# This script computes basic descriptive statistics and exploratory summaries
# across four populations: Turkana, Orang Asli, Tsimane, and Tarahumara.
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Complete case counts for each pain outcome × covariates
#      (Turkana and Orang Asli)
#   2. Density plots for key variables (industrialization index, BMI,
#      waist circumference, body fat, CRP, step counts, MVPA) —
#      both populations overlaid in a faceted ggplot
#   3. Descriptive summaries (five-number summary + NA count via
#      summarise_numeric(); sex breakdown by group)
#   4. Raw (unadjusted) pain prevalence tables for Turkana and Orang Asli
#      using summarise_binary_outcomes(), plus single-number prevalence
#      for Tsimane back pain and Tarahumara LBP
#
# Dependencies: loaded in 00b_setup.R
# Helper functions: defined in 00a_functions.R
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
# 1. COMPLETE CASE COUNTS
# ==============================================================================

# Complete cases for each pain outcome × covariates in Turkana and Orang Asli
count_complete_cases(dat_turkana)
count_complete_cases(dat_orangasli)

# ==============================================================================
# 2. DENSITY PLOTS FOR KEY VARIABLES
# ==============================================================================

# Variables to compare across both populations
density_vars <- c("industrialization_index", "bmi", "waist_circum_cm",
                  "body_fat_percentage", "crp_mg_l",
                  "step_count_average_daily", "mvpa_average_daily_min")

# Combine Turkana and Orang Asli, reshape to long format, facet by variable
bind_rows(
  dat_turkana  |> mutate(group = "Turkana"),
  dat_orangasli |> mutate(group = "Orang Asli")
) |>
  select(group, all_of(density_vars)) |>
  pivot_longer(-group, names_to = "variable", values_to = "value") |>
  drop_na(value) |>
  ggplot(aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.2, linewidth = 0.3) +
  facet_wrap(~ variable, scales = "free") +
  labs(x = NULL, y = "Density", fill = NULL, colour = NULL) +
  theme_bw_lbp() +
  theme(legend.position = "top")

# ==============================================================================
# 3. DESCRIPTIVE SUMMARIES
# ==============================================================================

# Five-number summaries + NA count for key numeric variables by group
summarise_numeric(
  Turkana      = dat_turkana,
  `Orang Asli` = dat_orangasli,
  vars = c("age_years", "industrialization_index", "bmi", "waist_circum_cm",
           "body_fat_percentage", "crp_mg_l",
           "step_count_average_daily", "mvpa_average_daily_min")
)

# Sex breakdown by group
dat_turkana |>
  count(sex) |>
  mutate(group = "Turkana", pct = n / sum(n) * 100)

dat_orangasli |>
  count(sex) |>
  mutate(group = "Orang Asli", pct = n / sum(n) * 100)

# ==============================================================================
# 4. RAW (UNADJUSTED) PAIN PREVALENCE
# ==============================================================================

# Prevalence tables for Turkana and Orang Asli via summarise_binary_outcomes(),
# plus single-number back pain prevalence for Tsimane and Tarahumara

# prevalence Turkana
summarise_binary_outcomes(
  dat_turkana,
  pattern = "pain.*(ny_0?1)$",
  label_fn = nice_labels
)

#   outcome               percent
# 1 joint pain 01            50  
# 2 knee pain 01             28  
# 3 low back pain n y 0 1    20.7
# 4 thoracic back pain 01    20.1
# 5 shoulder pain 01         15.6
# 6 foot pain 01             14.6
# 7 hip pain 01              13.8
# 8 elbow pain 01            11.9
# 9 neck pain 01             11.1
# 10 hand pain 01             7.5
# 11 widespread pain 01       6.7         

# prevalence Orang Asli
summarise_binary_outcomes(
  dat_orangasli,
  pattern = "pain.*(ny_0?1)$",
  label_fn = nice_labels
)

#   outcome                    percent
# 1 joint pain n y 0 1            64.2
# 2 low back pain n y 0 1         32.1
# 3 knee pain n y 0 1             28.6
# 4 shoulder pain n y 0 1         26  
# 5 thoracic back pain n y 0 1    20  
# 6 neck pain n y 0 1             19.4
# 7 hand pain n y 0 1             17.2
# 8 widespread pain n y 0 1       16.7
# 9 foot pain n y 0 1             15.9
# 10 elbow pain n y 0 1           14.8
# 11 hip pain n y 0 1             12.1          

# prevalence Tsimane
mean(dat_tsimane$back_pain, na.rm=TRUE) # 31.4% have back pain

# prevalence Tarahumara
mean(dat_tarahumara$low_back_pain, na.rm=TRUE) # 19.7% have low back pain


