
################################################################################
#
# DESCRIPTIVE ANALYSES — TARAHUMARA
# Low back pain prevalence (single outcome: low_back_pain, males only)
#
# This script computes descriptive prevalence statistics for low back pain
# among the Tarahumara. The sample comprises males only, so no sex
# stratification or sex comparisons are performed. Analyses are inline
# (not using the multi-group shared functions, which assume sex exists).
#
# The script is organized into the following sections:
#   0. Setup & data preparation
#   1. Age distribution and raw prevalence
#   2. Age-prevalence curve and age-adjusted prevalence (GAM)
#   3. Linearity test (GAM vs linear nested comparison)
#   4. Plots (age-adjusted prevalence + age-prevalence curve)
#
# Statistical approach:
#   - GAM with binomial(logit) link and penalized thin-plate spline
#     (k=6, REML, select=TRUE).
#   - Age-adjustment uses marginal standardization over the common age
#     support (99th percentile cap via common_age_support()).
#   - CIs use Wald intervals on the link scale, back-transformed via
#     inverse-logit (bounded, stable).
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

source("code/00b_setup.R")

set.seed(02138)


# --- Common age support (data-driven upper cap) ---
age_support <- common_age_support(dat_tarahumara, upper_quantile = 0.99, step = 1)
age_grid <- age_support$age_grid



# ==============================================================================
# 1. AGE DISTRIBUTION AND RAW PREVALENCE
# ==============================================================================

# age distribution
plot(density(dat_tarahumara$age_years))

# raw prevalence by age decade
dat_tarahumara |>
  group_by(age_bins = cut(age_years, breaks = c(10,20,30,40,50,60,70,80,90,100))) |>
  summarise(prev = mean(low_back_pain, na.rm=TRUE), n = n())



# ==============================================================================
# 2. AGE-PREVALENCE CURVE AND AGE-ADJUSTED PREVALENCE
# ==============================================================================

# Fit GAM (no sex term — males only)
m_gam <- gam(
  low_back_pain ~ s(age_years, k = 6, bs = "tp"),
  data   = dat_tarahumara,
  family = binomial(link = "logit"),
  method = "REML",
  select = TRUE
)

# Smooth p-value for the age effect 
summary(m_gam)$s.table

# Prediction grid
nd <- datagrid(model = m_gam, age_years = age_grid)

z <- qnorm(0.975)

# --- Age-prevalence curve data (Wald CIs on link scale, back-transformed) ---
avg_predictions(
  m_gam,
  newdata = nd,
  by      = "age_years",
  type    = "link"
) |>
  mutate(
    .eta      = estimate,
    estimate  = plogis(.eta) * 100,
    conf.low  = plogis(.eta - z * std.error) * 100,
    conf.high = plogis(.eta + z * std.error) * 100
  ) |>
  select(age_years, estimate, conf.low, conf.high) ->
prev_by_age

# --- Age-adjusted prevalence (marginal standardization over age grid, %) ---
avg_predictions(
  model   = m_gam,
  newdata = nd,
  type    = "response"
) |>
  mutate(
    estimate  = estimate  * 100,
    conf.low  = conf.low  * 100,
    conf.high = conf.high * 100
  ) ->
prev_age_adj

cat("Age-adjusted LBP prevalence (Tarahumara, males):\n")
print(prev_age_adj)


# ==============================================================================
# 4. PLOTS
# ==============================================================================

# Panel A: age-adjusted prevalence (single point + CI)
prev_age_adj |>
  ggplot(aes(x = "Male", y = estimate)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = NULL, y = "Age-Adjusted Prevalence", colour = NULL) +
  theme_bw_lbp() ->
p_adj

# Panel B: prevalence by age
prev_by_age |>
  ggplot(aes(x = age_years, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  geom_rug(data = dat_tarahumara |> filter(age_years >= age_support$bounds$lo,
                                           age_years <= age_support$bounds$hi),
           aes(x = age_years, y = 0), sides = "b",
           linewidth = 0.1, length = grid::unit(0.035, "npc"), inherit.aes = FALSE,
           position = position_jitter(width = 0.25, height = 0, seed = 123),
           alpha = 0.7, na.rm = TRUE) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(x = "Age (Years)", y = "Prevalence of Back Pain", colour = NULL, fill = NULL) +
  theme_bw_lbp() ->
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
       file = file.path(base_dir, "figures/tarahumara_all_pain.pdf"),
       width = 2.5, height = 3.4, units = "in")

