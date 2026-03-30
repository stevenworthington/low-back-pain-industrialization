
################################################################################
#
# DESCRIPTIVE ANALYSES — TSIMANE
# Back pain prevalence (single outcome: back_pain)
#
# This script computes descriptive prevalence statistics for back pain
# among the Tsimane, using the same shared functions as the multi-group
# analysis (02_descriptive_orangasli_turkana.R). Cluster-robust standard errors
# (by participant id) are used throughout because of repeated observations.
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Age distribution and raw prevalence
#   2. Age-adjusted prevalence (by sex, pooled, sex PR & difference)
#   3. Back pain age-prevalence curves by sex
#   4. Within-group sex-curve test (nested GAM comparison)
#   5. Plots (age-adjusted prevalence + age-prevalence curves)
#   6. Primary cause of pain (Tsimane-specific)
#
# Statistical approach:
#   - GAMs with binomial(logit) link and sex-specific penalized thin-plate
#     splines (k=6, REML, select=TRUE).
#   - Age-adjustment uses marginal standardization over the common age
#     support (99th percentile cap).
#   - Cluster-robust SEs (by participant id) via vcov = ~id passed to
#     marginaleffects functions.
#   - Sex PRs and differences via marginaleffects (avg_comparisons).
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

source("code/0b_setup.R")

set.seed(02138)


# Common age support (data-driven upper cap)
age_support <- common_age_support(dat_tsimane, upper_quantile = 0.99, step = 1)
age_ref  <- age_support$age_ref
age_grid <- age_support$age_grid



# ==============================================================================
# 1. AGE DISTRIBUTION AND RAW PREVALENCE
# ==============================================================================

# age distribution
plot(density(dat_tsimane$age_years))

# raw prevalence by age decade
dat_tsimane |>
  group_by(age_bins = cut(age_years, breaks = c(10,20,30,40,50,60,70,80))) |>
  summarise(prev = mean(back_pain, na.rm=TRUE), n = n())



# ==============================================================================
# 2. AGE-ADJUSTED PREVALENCE (BY SEX, POOLED, SEX PR & DIFFERENCE)
# ==============================================================================

# fit_predict_one() fits a sex-specific GAM per group and returns:
#   by_sex    – age-standardized prevalence by sex
#   pooled    – age- and sex-standardized prevalence (equal sex weights)
#   age_pvals – smooth p-values for the age effect by sex
#   sex_pr    – within-group Male/Female prevalence ratio (H0: PR = 1)
#   sex_diff  – within-group Male - Female prevalence difference
#
# Cluster-robust SEs are used via vcov = ~id.

fit_predict_one(
  outcome_var = "back_pain",
  df          = dat_tsimane,
  age_grid    = age_ref,
  pool        = "equal",
  vcov        = ~ id
  ) ->
tsimane_prev

# Age-adjusted prevalence by sex (%)
tsimane_prev$by_sex

# Age- and sex-standardized overall prevalence (%)
tsimane_prev$pooled

# Smooth p-values for age effect
tsimane_prev$age_pvals

# Sex prevalence ratio (Male / Female)
tsimane_prev$sex_pr

# Sex prevalence difference (Male - Female, %)
tsimane_prev$sex_diff



# ==============================================================================
# 3. BACK PAIN AGE-PREVALENCE CURVES BY SEX
# ==============================================================================

# Uses prevalence_by_age_x_sex_curve() for GAM-predicted prevalence (%)
# by age × sex, with Wald CIs on the link scale back-transformed via
# inverse-logit. Cluster-robust SEs via vcov = ~id.

prevalence_by_age_x_sex_curve(
  outcomes       = "back_pain",
  df             = dat_tsimane,
  k              = 6,
  bs             = "tp",
  upper_quantile = 0.99,
  age_step       = 1,
  vcov           = ~ id
  ) ->
curve_df



# ==============================================================================
# 4. WITHIN-GROUP SEX-CURVE TEST (NESTED GAM COMPARISON)
# ==============================================================================

# Tests whether the age-prevalence curve shape differs by sex.
#   m0: common age smooth (sex as main effect only)
#   m1: sex-specific age smooths
# Uses ML (not REML) for a valid LRT.

m0 <- gam(
  back_pain ~ sex + s(age_years, k = 6, bs = "tp"),
  data   = dat_tsimane,
  family = binomial(link = "logit"),
  method = "ML",
  select = FALSE
)

m1 <- gam(
  back_pain ~ sex + s(age_years, by = sex, k = 6, bs = "tp"),
  data   = dat_tsimane,
  family = binomial(link = "logit"),
  method = "ML",
  select = FALSE
)

anova(m0, m1, test = "Chisq")
# p = 0.92



# ==============================================================================
# 5. PLOTS
# ==============================================================================

# Panel A: sex-specific age-adjusted prevalence
tsimane_prev$by_sex |>
  ggplot(aes(x = sex, y = estimate, colour = sex)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = NULL, y = "Age-Adjusted Prevalence", colour = NULL) +
  theme_bw_lbp() +
  theme(legend.position = "none") ->
p_adj

# Panel B: prevalence by sex × age 
curve_df |>
  ggplot(aes(x = age_years, y = estimate, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  geom_rug(data = dat_tsimane |> filter(age_years >= age_support$bounds$lo,
                                        age_years <= age_support$bounds$hi),
           aes(x = age_years, colour = sex, y = 0), sides = "b",
           linewidth = 0.1, length = grid::unit(0.035, "npc"), inherit.aes = FALSE,
           position = position_jitter(width = 0.25, height = 0, seed = 123),
           alpha = 0.7, na.rm = TRUE) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = "Age (Years)", y = "Prevalence of Back Pain", colour = NULL, fill = NULL) +
  theme_bw_lbp() +
  theme(legend.position = "top") ->
p_curve

# Combined figure
p_adj_mod   <- p_adj + guides(colour = "none", fill = "none",
                               linetype = "none", shape = "none")
p_curve_mod <- p_curve + theme(strip.text = element_blank())

(p_adj_mod / plot_spacer() / p_curve_mod) +
  plot_layout(ncol = 1, guides = "collect", heights = c(1, -0.05, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.margin     = margin(b = -8),
    legend.box.margin = margin(b = 0),
    plot.tag          = element_text(size = 10, face = "plain"),
    plot.tag.position = c(0.012, 0.985),
    plot.margin       = margin(2, 2, 2, 2),
    axis.title.y = element_text(margin = margin(r = 4))
  ) ->
p_combined

p_combined <- p_combined + theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"))

ggsave(plot = p_combined,
       file = file.path(base_dir, "figures/tsimane_all_pain.pdf"),
       width = 2.5, height = 3.5, units = "in")



# ==============================================================================
# 6. PRIMARY CAUSE OF PAIN (TSIMANE-SPECIFIC)
# ==============================================================================

# Among people who reported back pain. primary_cause is only non-NA for
# people with back pain; NA = not asked / no back pain / missing → exclude.

dat_tsimane |>
  filter(!is.na(primary_cause)) |>
  mutate(primary_cause = str_to_sentence(primary_cause),
         primary_cause = fct_infreq(primary_cause)) ->
cause_df

cause_df |>
  count(primary_cause) |>
  mutate(prop = n / sum(n)) ->
cause_prop

# Proportions bar chart 
cause_prop |>
ggplot(aes(x = prop, y = primary_cause)) +
  geom_col(fill = "grey60", width = 0.9) +
  geom_text(aes(label = paste0(n, " (", scales::percent(prop, accuracy = 1), ")")),
            hjust = -0.1, size = 2) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.23))) +
  labs(x = "Proportion of Reported Primary Causes",
       y = "Self-Reported\nPrimary Cause of Back Pain") +
  theme_bw_lbp() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(r = 4))) ->
cause_prop_plot

ggsave(plot = cause_prop_plot,
       file = file.path(base_dir, "figures/tsimane_cause_prop_barplot.pdf"),
       width = 3.5, height = 2, units = "in")

# Counts bar chart
cause_prop |>
ggplot(aes(x = n, y = primary_cause)) +
  geom_col(fill = "grey60", width = 0.9) +
  geom_text(aes(label = paste0(n, " (", scales::percent(prop, accuracy = 1), ")")),
            hjust = -0.1, size = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.23))) +
  labs(x = "Number of Respondents",
       y = "Self-Reported\nPrimary Cause of Back Pain") +
  theme_bw_lbp() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(r = 4))) ->
cause_count_plot

ggsave(plot = cause_count_plot,
       file = file.path(base_dir, "figures/tsimane_cause_count_barplot.pdf"),
       width = 3.5, height = 2, units = "in")

