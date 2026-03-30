
################################################################################
#
# DESCRIPTIVE ANALYSES
# Multi-group pain prevalence: Turkana & Orang Asli
#
# This script computes descriptive prevalence statistics for musculoskeletal
# pain outcomes across two populations (Turkana and Orang Asli), broken down
# by sex and age.
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Age distributions by sex and group
#   2. Age-adjusted prevalence by pain location, sex, and group
#      (dot + CI plot, pooled estimates, within/between-group PRs & diffs)
#   3. LBP age-prevalence curves by sex and group (ribbon plot)
#   4. Combined two-panel figure (patchwork)
#   5. Within-group sex-curve test for LBP (nested GAM comparison)
#   6. Age effect tests, monotonicity, plateau ages, and group × age interactions
#   7. Supplemental: age-prevalence curves for all pain locations
#
# Statistical approach:
#   - GAMs with binomial(logit) link and sex-specific penalized thin-plate
#     splines (k=6, REML, select=TRUE) are the workhorse throughout.
#   - Age-adjustment uses marginal standardization over observed ages within
#     the common age support across groups (avoids extrapolation).
#   - CIs for age-specific curves use Wald intervals on the link scale,
#     back-transformed via inverse-logit (bounded, stable).
#   - Within-group sex PRs and differences use marginaleffects
#     (avg_comparisons with comparison="ratio" and "difference").
#   - Between-group PRs use the delta method on the log scale from
#     independent pooled estimates.
#
# Dependencies: loaded in 0b_setup.R
# Helper functions: defined in 0a_functions.R
#
################################################################################


# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================

base_dir <- "."
setwd(base_dir)

# Source setup file (loads packages, data, custom functions, themes)
source("code/0b_setup.R")

set.seed(02138)


# --- Subset to complete cases for modeling ---
# Keep only rows with non-missing outcome, sex, age, and group
dat <- dat |>
  filter(!is.na(low_back_pain_ny_01), !is.na(sex), !is.na(age_years), !is.na(group))


# --- Outcome definitions ---
# All 11 binary pain outcomes and their display labels
outcomes <- c(
  "widespread_pain_ny_01", "joint_pain_ny_01", "neck_pain_ny_01",
  "shoulder_pain_ny_01", "thoracic_back_pain_ny_01",
  "low_back_pain_ny_01", "hip_pain_ny_01", "elbow_pain_ny_01",
  "hand_pain_ny_01", "knee_pain_ny_01", "foot_pain_ny_01"
)

levs <- c(
  "Widespread", "Any Joint", "Neck", "Shoulder", "Thoracic Back",
  "Low Back", "Hip", "Elbow", "Hand", "Knee", "Foot"
)

lab_map <- setNames(levs, outcomes)


# --- common age support across groups ---
# uses common_age_support() from 0a_functions.R to compute the overlapping
# age range. The upper bound is data-driven (99th percentile of pooled ages,
# capped at the smallest group maximum). age_ref contains observed ages for
# standardization; age_grid is a regular 1-year sequence for curve plotting.
age_support <- common_age_support(dat, upper_quantile = 0.99, step = 1)
age_ref  <- age_support$age_ref
age_grid <- age_support$age_grid



# ==============================================================================
# 1. AGE DISTRIBUTIONS BY SEX AND GROUP
# ==============================================================================

# --- ECDF by sex, faceted by group ---
# Empirical CDFs reveal full-distribution differences in age sampling
ggplot(dat, aes(age_years, colour = sex)) +
  stat_ecdf(linewidth = 0.7) +
  facet_wrap(~ group) +
  labs(x = "Age (years)", y = "ECDF") +
  theme_bw()

# --- Density by sex, faceted by group ---
ggplot(dat, aes(age_years, fill = sex, colour = sex)) +
  geom_density(alpha = 0.2, linewidth = 0.3) +
  facet_wrap(~ group) +
  labs(x = "Age (years)", y = "Density") +
  theme_bw()

# --- Summary table: age by group × sex ---
dat |>
  group_by(group, sex) |>
  summarise(
    n      = n(),
    mean   = mean(age_years),
    sd     = sd(age_years),
    median = median(age_years),
    IQR    = IQR(age_years),
    min    = min(age_years),
    max    = max(age_years),
    .groups = "drop"
  ) ->
sumtab

# raw prevalence by age decade
dat |>
  group_by(group, 
           age_bins = cut(age_years, breaks = c(10,20,30,40,50,60,70,80,90,100))) |>
  summarise(prev = mean(low_back_pain_ny_01, na.rm=TRUE), n = n())



# ==============================================================================
# 2. AGE-ADJUSTED PREVALENCE BY PAIN LOCATION, SEX, AND GROUP
# ==============================================================================

# uses fit_predict_all() from 0a_functions.R, which wraps fit_predict_one()
# across all outcomes. For each outcome × group, it fits a sex-specific GAM
# and computes several estimands via marginal standardization over the common
# age support, with marginaleffects for ratio and difference inference.
#
# Returns a named list of seven tibbles (each with a `location` column):
#   by_sex            – age-standardized prevalence by sex (primary estimand)
#   pooled            – age- and sex-standardized prevalence per group (pooled across
#                       sex with equal weights because pool="equal")
#   age_pvals         – smooth p-values for the age effect per group × sex
#   sex_pr            – within-group Male/Female prevalence ratios (H0: PR = 1)
#   sex_diff          – within-group Male - Female prevalence differences
#   group_contrast_pr – between-group prevalence ratios with CI, p-value, and % change

fit_predict_all(  
  outcomes    = outcomes,
  lab_map     = lab_map,
  df          = dat,
  age_grid    = age_ref,
  pool        = "equal",
  group_order = c("Orang Asli", "Turkana")
  ) ->
all_prev

# compute overall mean prevalence per location to order the y-axis
all_prev$by_sex |>
  group_by(location) |>
  summarise(m = mean(estimate, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(m)) |>
  pull(location) ->
location_order

# dot + CI plot: prevalence by location, sex, group 
all_prev$by_sex |>
  ggplot(aes(x = estimate, y = location, colour = sex)) +
  geom_point(position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0,
                position = position_dodge(width = 0.6)) +
  facet_wrap(~ group, ncol = 2) +
  scale_x_continuous(labels = scales::label_percent(scale = 1)) +
  scale_y_discrete(limits = location_order) +
  labs(x = "Age-Adjusted Prevalence", y = "Pain Location", colour = NULL) +
  theme_bw_lbp() +
  theme(legend.position.inside = c(0.4, 0.9)) ->
p1_age_standardized_prevalence_by_sex



# ==============================================================================
# 3. LBP AGE-PREVALENCE CURVES BY SEX AND GROUP
# ==============================================================================

# fit per-group GAMs and extract age-specific prevalence predictions.
# CIs are computed on the link (logit) scale and back-transformed via
# inverse-logit for bounded [0, 1] intervals (Wald approach).

prevalence_by_age_x_sex_curve(
  outcomes       = "low_back_pain_ny_01",
  df             = dat,
  lab_map        = lab_map,
  k              = 6,
  bs             = "tp",
  upper_quantile = 0.99
  ) ->
prev_lbp_age

# rug data for age distribution overlay 
rug_df <- dat |> filter(age_years >= lo, age_years <= hi)

# age-prevalence curve plot 
prev_lbp_age |>
  ggplot(aes(x = age_years, y = estimate, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  geom_rug(data = rug_df, aes(x = age_years, colour = sex, y = 0), sides = "b",
           linewidth = 0.1, length = grid::unit(0.035, "npc"), inherit.aes = FALSE,
           position = position_jitter(width = 0.25, height = 0, seed = 123),
           alpha = 0.7, na.rm = TRUE) +
  facet_wrap(~ group) +
  coord_cartesian(xlim = c(lo, hi)) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = "Age (Years)", y = "Prevalence of LBP", colour = NULL, fill = NULL) +
  theme_bw_lbp() +
  theme(legend.position.inside = c(0.6, 0.9)) ->
p2_age_prevalence_curve



# ==============================================================================
# 4. COMBINED FIGURE (PATCHWORK)
# ==============================================================================

# Panel A: prevalence by location (dot + CI)
# Panel B: LBP age-prevalence curves
# Shared legend at top; strip labels removed from bottom panel (redundant)

p1_mod <- p1_age_standardized_prevalence_by_sex + 
              guides(colour = "none", fill = "none",
                     linetype = "none", shape = "none")
p2_mod <- p2_age_prevalence_curve + 
              theme(strip.text = element_blank())

final <- (p1_mod + plot_spacer() + p2_mod) +
  plot_layout(ncol = 1, guides = "collect", heights = c(1.5, -0.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.margin     = margin(b = -14),
    legend.box.margin = margin(b = 0),
    plot.tag          = element_text(size = 10, face = "plain"),
    plot.tag.position = c(0.012, 0.985),
    plot.margin       = margin(2, 2, 2, 2)
  )

final <- final + theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"))

ggsave(plot = final,
       file = file.path(base_dir, "figures/combined_lbp_plots.pdf"),
       width = 4.5, height = 3.5, units = "in")



# ==============================================================================
# 5. WITHIN-GROUP SEX-CURVE TEST FOR LBP (NESTED GAM COMPARISON)
# ==============================================================================

# Tests whether LBP age-prevalence curves differ by sex within each group
# using a nested GAM comparison (likelihood ratio test).
#   m0: common age smooth for both sexes (sex as main effect only)
#   m1: sex-specific age smooths (allows shape to vary by sex)
# Uses ML (not REML) for a valid LRT.

groups_list <- split(dat, dat$group)

sex_curve_tests <- lapply(groups_list, function(df_g) {

  grp_name <- as.character(df_g$group[1])

  m0 <- gam(
    low_back_pain_ny_01 ~ sex + s(age_years, k = 6, bs = "tp"),
    data = df_g, family = binomial(link = "logit"),
    method = "ML", select = FALSE
  )
  
  m1 <- gam(
    low_back_pain_ny_01 ~ sex + s(age_years, k = 6, bs = "tp") +
      s(age_years, by = sex, k = 6, bs = "tp"),
    data = df_g, family = binomial(link = "logit"),
    method = "ML", select = FALSE
  )

  list(
    group      = grp_name,
    curve_test = anova(m0, m1, test = "Chisq")
  )
})

names(sex_curve_tests) <- names(groups_list)

# Print LRT results
for (nm in names(sex_curve_tests)) {
  cat("\n===", nm, "===\n")
  cat("Age-curve sex interaction (LRT):\n")
  print(sex_curve_tests[[nm]]$curve_test)
}



# ==============================================================================
# 6. AGE EFFECT TESTS, MONOTONICITY, PLATEAU AGES, AND GROUP × AGE INTERACTIONS
# ==============================================================================

# Three complementary summaries for LBP:
#   a) Smooth p-values: is the age-prevalence relationship non-flat?
#      (already computed in Section 2 via fit_predict_one → all_prev$age_pvals)
#   b) Monotonicity: is the derivative significantly > 0 across the age grid?
#      (from derivative_plateau_slopes(), which also provides plateau ages in §7)
#   c) Group × age interaction: do age-prevalence curve shapes differ across
#      groups within each sex? (nested GAM LRT via test_group_x_age())

# a) Smooth p-values from Section 2
#    Filter to LBP; for other outcomes, all_prev$age_pvals has the full set.
all_prev$age_pvals |>
  filter(location == "Low Back") ->
lbp_age_pvals

cat("Age effect p-values by group × sex (LBP):\n")
print(lbp_age_pvals)
# Use p_age for "prevalence varied with age" claims per population × sex.

# b) Monotonicity and plateau ages via derivative_plateau_slopes()
dat |>
  split(~group) |>
  purrr::map_dfr(function(df_g) {
    purrr::map_dfr(c("Female", "Male"), function(sx) {
      res <- derivative_plateau_slopes(df_g, outcome = "low_back_pain_ny_01", sex_level = sx)
      tibble(group = as.character(df_g$group[1]),
             sex = res$sex, increasing = res$increasing,
             frac_positive = res$frac_positive, plateau_age = res$plateau_age)
    })
  }) ->
deriv_results

cat("\nMonotonicity and plateau ages (LBP):\n")
print(deriv_results)
# increasing:    TRUE if derivative CI lower bound > 0 across the entire age grid
# frac_positive: fraction of the age grid where the derivative is significantly > 0
# plateau_age:   upper bound of the longest contiguous positive-derivative run
#                (the age at which prevalence stops significantly increasing)

# c) Group × age interaction (nested GAM LRT per sex)
test_group_x_age(
  data    = dat,
  outcome = "low_back_pain_ny_01",
  k       = 6,
  bs      = "tp"
  ) ->
group_x_age
  
cat("\nGroup × age interaction p-values by sex (LBP):\n")
print(group_x_age)
# Use p_group_x_age for group × age interaction within females and males
# (e.g., "steeper among Orang Asli females, p = ...").



# ==============================================================================
# 7. SUPPLEMENTAL: AGE-PREVALENCE CURVES FOR ALL PAIN LOCATIONS
# ==============================================================================

# build a faceted grid of age-prevalence curves for all pain outcomes
# (excluding LBP which is shown in the main figure).

outcomes_supp <- setdiff(outcomes, "low_back_pain_ny_01")

pred_df <- prevalence_by_age_x_sex_curve(
  outcomes       = outcomes_supp,
  df             = dat,
  lab_map        = lab_map,
  k              = 6,
  bs             = "tp",
  upper_quantile = 0.99,
  age_step       = 1
)

# supplemental figure: prevalence by age, sex, location, group
pred_df |>
  mutate(label = factor(label, levels = setdiff(location_order, "Low Back"))) |>
  ggplot(aes(age_years, estimate, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_grid(label ~ group, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = "Age (Years)", y = "Pain Prevalence by Location",
       colour = NULL, fill = NULL) +
  theme_bw_lbp() +
  theme(
    legend.position = "top",
    legend.box.margin = margin(b = -12),
    legend.margin = margin(t = 0, b = 0),
    strip.placement = "outside",
    strip.text.y.left = element_text(margin = margin(l = 5, r = 1)),
    strip.switch.pad.grid = unit(1, "pt")
  ) ->
prob_all_pain

ggsave(plot = prob_all_pain,
       file = file.path(base_dir, "figures/prob_all_pain.pdf"),
       width = 4.7, height = 8, units = "in")


