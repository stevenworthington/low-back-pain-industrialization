################################################################################
#
# Pr(LBP) vs RISK FACTORS -- VISUALIZATION
# Low back pain causal analysis pipeline
#
# This script visualizes the causal effects of six risk factors on the
# probability of low back pain (LBP). The six risk-factor exposures are:
# BMI, waist circumference, body fat %, C-reactive protein (CRP, log scale),
# mean daily step counts, and moderate-to-vigorous physical activity (MVPA).
#
# For each risk factor, the script:
#   - Combines per-group AERF data (Turkana, Orang Asli) with the between-
#     group difference AERF (from clarify) into a single tibble. Per-group
#     values are on the probability scale and are multiplied by 100 to yield
#     percent. Difference data uses different column names (Estimate, 2.5 %,
#     97.5 %) so requires renaming and vctrs::vec_data() to strip class
#     attributes.
#   - Repeats the same combination for AMEF data (also x100 for percent).
#   - Calls plot_AERF_AMEF() with outcome = quote(Pr(LBP)) and risk-factor-
#     specific axis labels and limits.
#
# After individual panels are built, the script assembles them into composite
# figures:
#   - A 3x2 grid of per-group panels with shared legend, saved to
#     figures/pain_plots_groups.pdf.
#   - A 3x2 grid of difference panels, saved to figures/pain_plots_diffs.pdf.
#
# Sections:
#   0. Setup -- source helpers, set seed, load saved model .RData files
#   1. Pr(LBP) vs BMI
#   2. Pr(LBP) vs waist circumference
#   3. Pr(LBP) vs body fat %
#   4. Pr(LBP) vs CRP (log-scale x-axis via transform_expo="log")
#   5. Pr(LBP) vs daily step counts
#   6. Pr(LBP) vs MVPA
#   7. Combine per-group panels into 3x2 grid figure
#   8. Combine between-group difference panels into 3x2 grid figure
#
# Dependencies:
#   00b_setup.R (packages, helpers, data); saved model objects in models/
#   Packages: tidyverse, patchwork, cowplot, vctrs, ggplot2
#   Helper functions: plot_AERF_AMEF(), get_legend() (cowplot)
#
# Outputs:
#   figures/pain_plots_groups.pdf  -- 5 x 8 in composite (per-group panels)
#   figures/pain_plots_diffs.pdf   -- 5 x 8 in composite (difference panels)
#
################################################################################


# ==============================================================================
# 0. SETUP
# ==============================================================================

base_dir <- "."
setwd(base_dir)

# source setup file
source("code/00b_setup.R")

set.seed(02138)


################################################################################
##### load data
################################################################################

# get filenames for every .RData in /data
models_filenames <- list.files(path = "models", pattern = ".RData", 
                               full.names = TRUE, ignore.case = TRUE)

# load every .RData in /data into the global environment
walk(models_filenames, ~ load(.x, envir = .GlobalEnv))


# ==============================================================================
# 1. Pr(LBP) vs BMI
# ==============================================================================

# Combine per-group and difference AERF data. Per-group estimates are on the
# probability scale, so multiply by 100 to express as percentages.
bind_rows(
  pain_risk_bmi_turkana$AERF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_bmi_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_bmi_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_bmi_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_bmi_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_bmi_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_bmi_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_bmi_AMEF
  
# Generate paired AERF + AMEF panels for Pr(LBP) vs BMI
pain_risk_bmi_plots <- plot_AERF_AMEF(
  dat_AERF   = pain_risk_bmi_AERF, 
  dat_AMEF   = pain_risk_bmi_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = quote("pp/"~kg/m^2),
  legend_pos = c(0.25, 0.85),
  exposure   = "bmi",
  xlab       = expression("Body Mass Index (kg/m"^2*")"), 
  xlim       = c(13.6, 42.3)
  )

# save
# ggsave(plot=pain_risk_bmi_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_bmi_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_bmi_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_bmi_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 2. Pr(LBP) vs WAIST CIRCUMFERENCE
# ==============================================================================

# Combine per-group and difference AERF data (x100 for percent)
bind_rows(
  pain_risk_waist_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_waist_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_waist_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_waist_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_waist_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_waist_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_waist_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_waist_AMEF

# Generate paired AERF + AMEF panels for Pr(LBP) vs waist circumference
pain_risk_waist_plots <- plot_AERF_AMEF(
  dat_AERF   = pain_risk_waist_AERF, 
  dat_AMEF   = pain_risk_waist_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = "pp/cm",
  legend_pos = c(0.25, 0.85),
  exposure   = "waist_circum_cm",
  xlab       = "Waist Circumference (cm)", 
  xlim       = c(53.8, 116.8)
  )

# save
# ggsave(plot=pain_risk_waist_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_waist_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_waist_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_waist_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 3. Pr(LBP) vs BODY FAT %
# ==============================================================================

# Combine per-group and difference AERF data (x100 for percent)
bind_rows(
  pain_risk_fat_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_fat_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_fat_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_fat_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_fat_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_fat_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_fat_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_fat_AMEF

# Generate paired AERF + AMEF panels for Pr(LBP) vs body fat %
pain_risk_fat_plots <- plot_AERF_AMEF(
  dat_AERF   = pain_risk_fat_AERF, 
  dat_AMEF   = pain_risk_fat_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = "pp/%-pt",
  legend_pos = c(0.25, 0.85),
  exposure   = "body_fat_percentage",
  xlab       = "Body Fat (%)", 
  xlim       = c(4.9, 49.1)
 )

# save
# ggsave(plot=pain_risk_fat_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_fat_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_fat_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_fat_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 4. Pr(LBP) vs CRP (log-scale x-axis)
# ==============================================================================

# Combine per-group and difference AERF data (x100 for percent)
bind_rows(
  pain_risk_crp_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_crp_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_crp_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_crp_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_crp_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_crp_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_crp_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_crp_AMEF

# Generate paired AERF + AMEF panels; transform_expo="log" applies log x-axis
pain_risk_crp_plots <- plot_AERF_AMEF(
  dat_AERF       = pain_risk_crp_AERF, 
  dat_AMEF       = pain_risk_crp_AMEF,
  outcome        = quote(Pr(LBP)),
  units_AERF     = NULL,
  units_AMEF     = quote("pp/"~mg/L),
  legend_pos     = c(0.75, 0.85),
  transform_expo = "log",
  exposure       = "crp_mg_l",
  xlab           = "CRP (mg/L, log scale)", 
  xlim           = c(0.095, 57)
 )

# save
# ggsave(plot=pain_risk_crp_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_crp_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_crp_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_crp_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 5. Pr(LBP) vs DAILY STEP COUNTS
# ==============================================================================

# Combine per-group and difference AERF data (x100 for percent)
bind_rows(
  pain_risk_steps_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_steps_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_steps_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_steps_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_steps_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_steps_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_steps_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_steps_AMEF

# Generate paired AERF + AMEF panels for Pr(LBP) vs daily steps
pain_risk_steps_plots <- plot_AERF_AMEF(
  dat_AERF = pain_risk_steps_AERF, 
  dat_AMEF = pain_risk_steps_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = "pp/step",
  legend_pos = c(0.75, 0.85),
  exposure   = "step_count_average_daily",
  xlab       = "Mean Daily Steps", 
  xlim       = c(1900, 32600)
  )

# save
# ggsave(plot=pain_risk_steps_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_steps_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_steps_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_steps_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 6. Pr(LBP) vs MVPA
# ==============================================================================

# Combine per-group and difference AERF data (x100 for percent)
bind_rows(
  pain_risk_mvpa_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_mvpa_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_mvpa_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_mvpa_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_risk_mvpa_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_mvpa_orangasli$AMEF |>
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_risk_mvpa_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_risk_mvpa_AMEF

# Generate paired AERF + AMEF panels for Pr(LBP) vs MVPA
pain_risk_mvpa_plots <- plot_AERF_AMEF(
  dat_AERF = pain_risk_mvpa_AERF, 
  dat_AMEF = pain_risk_mvpa_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = "pp/min/day",
  legend_pos = c(0.75, 0.85),
  exposure   = "mvpa_average_daily_min",
  xlab       = "MVPA (min/day)", 
  xlim       = c(14, 405)
 )

# save
# ggsave(plot=pain_risk_mvpa_plots$group_plot,
#        file=file.path(base_dir, "figures/pain_mvpa_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=pain_risk_mvpa_plots$diff_plot,
#        file=file.path(base_dir, "figures/pain_mvpa_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 7. COMBINE PER-GROUP PANELS INTO 3x2 GRID
# ==============================================================================

# Helper: strip individual legends and tighten outer margins on each patchwork
# composite so panels sit close together in the final cowplot grid.
tighten_pair <- function(
    p,
    outer = margin(t=1, r=-3, b=1, l=-3, unit = "pt")) { # of the WHOLE pair
  (p + plot_annotation(theme = theme(plot.margin = outer))) &
       theme(legend.position = "none")
}

# Apply tighten_pair() to each risk factor's per-group panel
p1 <- tighten_pair(pain_risk_bmi_plots$group_plot)
p2 <- tighten_pair(pain_risk_waist_plots$group_plot)
p3 <- tighten_pair(pain_risk_fat_plots$group_plot)
p4 <- tighten_pair(pain_risk_crp_plots$group_plot)
p5 <- tighten_pair(pain_risk_steps_plots$group_plot)
p6 <- tighten_pair(pain_risk_mvpa_plots$group_plot)

# Extract a shared horizontal legend from the first panel
legend_plot <- p1 + guides(color = "none") +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(b = 2, t = 2))

leg <- get_legend(legend_plot)

# Assemble 3x2 grid (A-F labels) with shared legend on top
grid <- plot_grid(
  p1, p2,
  p3, p4,
  p5, p6,
  ncol = 2, align = "hv", axis = "tblr",
  labels = LETTERS[1:6], label_size = 10, label_fontface = "plain",
  label_x = 0.03, label_y = 0.985
)

final <- plot_grid(leg, grid, ncol = 1, rel_heights = c(0.03, 1))
final <- final + theme(plot.margin = margin(t=4,r=6,b=4,l=4, unit="pt"))

# Save per-group composite figure
ggsave(plot=final,
       file=file.path(base_dir, "figures/pain_plots_groups.pdf"),
       width = 5, height = 8, units = "in")


# ==============================================================================
# 8. COMBINE BETWEEN-GROUP DIFFERENCE PANELS INTO 3x2 GRID
# ==============================================================================

# Re-define tighten_pair() (identical to Section 7; redefined here for
# self-contained execution of this section)
tighten_pair <- function(
    p,
    outer = margin(t=1, r=-3, b=1, l=-3, unit = "pt")) { # of the WHOLE pair
  (p + plot_annotation(theme = theme(plot.margin = outer))) &
    theme(legend.position = "none")
}

# Apply tighten_pair() to each risk factor's difference panel
p1 <- tighten_pair(pain_risk_bmi_plots$diff_plot)
p2 <- tighten_pair(pain_risk_waist_plots$diff_plot)
p3 <- tighten_pair(pain_risk_fat_plots$diff_plot)
p4 <- tighten_pair(pain_risk_crp_plots$diff_plot)
p5 <- tighten_pair(pain_risk_steps_plots$diff_plot)
p6 <- tighten_pair(pain_risk_mvpa_plots$diff_plot)

# Assemble 3x2 grid (A-F labels); no shared legend needed for difference panels
grid <- plot_grid(
  p1, p2,
  p3, p4,
  p5, p6,
  ncol = 2, align = "hv", axis = "tblr",
  labels = LETTERS[1:6], label_size = 10, label_fontface = "plain",
  label_x = 0.03, label_y = 0.985
)

final <- grid + theme(plot.margin = margin(t=4,r=6,b=4,l=4, unit="pt"))

# Save between-group difference composite figure
ggsave(plot=final,
       file=file.path(base_dir, "figures/pain_plots_diffs.pdf"),
       width = 5, height = 8, units = "in")

