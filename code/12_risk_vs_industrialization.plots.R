################################################################################
#
# RISK FACTOR vs INDUSTRIALIZATION -- VISUALIZATION
# Low back pain causal analysis pipeline
#
# This script visualizes the causal effects of industrialization on six
# cardiometabolic and physical-activity risk factors: BMI, waist circumference,
# body fat percentage, C-reactive protein (CRP), mean daily step counts, and
# moderate-to-vigorous physical activity (MVPA).
#
# For each risk factor, the script:
#   - Combines per-group average estimated response function (AERF) data
#     (Turkana, Orang Asli) with the between-group difference AERF (produced
#     by the clarify package) into a single tibble. The difference data uses
#     different column names (Estimate, 2.5 %, 97.5 %) so requires renaming
#     and vctrs::vec_data() to strip class attributes.
#   - Repeats the same combination for average marginal effect function (AMEF)
#     data. CRP AMEF requires special back-transformation: per-group slopes
#     are on the log scale, so (exp(x)-1)*100 converts to percent change;
#     the difference is already on the percent-change scale from
#     create_pct_diff().
#   - Calls plot_AERF_AMEF() to produce paired AERF + AMEF panels for both
#     the per-group and between-group difference views.
#
# After producing individual panels, the script assembles them into composite
# figures:
#   - A 3x2 grid of per-group panels with a shared legend extracted via
#     cowplot::get_legend(), saved to figures/risk_plots_groups.pdf.
#   - A 3x2 grid of difference panels, saved to figures/risk_plots_diffs.pdf.
#   - Uses a tighten_pair() helper to strip legends and tighten margins on
#     each patchwork composite before grid assembly.
#
# Sections:
#   0. Setup -- source helpers, set seed, load saved model .RData files
#   1. BMI vs industrialization
#   2. Waist circumference vs industrialization
#   3. Body fat % vs industrialization
#   4. CRP vs industrialization (log-scale back-transformation for AMEF)
#   5. Daily step counts vs industrialization
#   6. MVPA vs industrialization
#   7. Combine per-group panels into 3x2 grid figure
#   8. Combine between-group difference panels into 3x2 grid figure
#
# Dependencies:
#   00b_setup.R (packages, helpers, data); saved model objects in models/
#   Packages: tidyverse, patchwork, cowplot, vctrs, ggplot2
#   Helper functions: plot_AERF_AMEF(), get_legend() (cowplot)
#
# Outputs:
#   figures/risk_plots_groups.pdf  -- 5 x 8 in composite (per-group panels)
#   figures/risk_plots_diffs.pdf   -- 5 x 8 in composite (difference panels)
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
# 1. BMI vs INDUSTRIALIZATION
# ==============================================================================

# Combine Turkana, Orang Asli, and between-group difference AERF into one tibble.
# Difference data from clarify uses different column names, so rename and strip
# class attributes with vctrs::vec_data() for compatibility with bind_rows().
bind_rows(
  risk_bmi_turkana$AERF |> mutate(group = "Turkana"),
  risk_bmi_orangasli$AERF |> mutate(group = "Orang Asli"),
  risk_bmi_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_bmi_AERF

# Combine AMEF data (same rename + strip pattern as AERF above)
bind_rows(
  risk_bmi_turkana$AMEF |> mutate(group = "Turkana"),
  risk_bmi_orangasli$AMEF |> mutate(group = "Orang Asli"),
  risk_bmi_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_bmi_AMEF
  
# Generate paired AERF + AMEF panels (per-group and difference views)
risk_bmi_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_bmi_AERF, 
  dat_AMEF   = risk_bmi_AMEF,
  outcome    = "BMI",
  units_AERF = quote(kg/m^2),
  units_AMEF = quote(kg/m^2),
  legend_pos = c(0.25, 0.85)
  )

# save
# ggsave(plot=risk_bmi_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_bmi_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_bmi_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_bmi_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 2. WAIST CIRCUMFERENCE vs INDUSTRIALIZATION
# ==============================================================================

# Combine per-group and difference AERF data
bind_rows(
  risk_waist_turkana$AERF |> mutate(group = "Turkana"),
  risk_waist_orangasli$AERF |> mutate(group = "Orang Asli"),
  risk_waist_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_waist_AERF

# Combine per-group and difference AMEF data
bind_rows(
  risk_waist_turkana$AMEF |> mutate(group = "Turkana"),
  risk_waist_orangasli$AMEF |> mutate(group = "Orang Asli"),
  risk_waist_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_waist_AMEF

# Generate paired AERF + AMEF panels
risk_waist_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_waist_AERF, 
  dat_AMEF   = risk_waist_AMEF,
  outcome    = "Waist Circ.",
  units_AERF = "cm",
  units_AMEF = "cm",
  legend_pos = c(0.25, 0.85)
  )

# save
# ggsave(plot=risk_waist_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_waist_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_waist_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_waist_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 3. BODY FAT % vs INDUSTRIALIZATION
# ==============================================================================

# Combine per-group and difference AERF data
bind_rows(
  risk_fat_turkana$AERF |> mutate(group = "Turkana"),
  risk_fat_orangasli$AERF |> mutate(group = "Orang Asli"),
  risk_fat_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_fat_AERF

# Combine per-group and difference AMEF data
bind_rows(
  risk_fat_turkana$AMEF |> mutate(group = "Turkana"),
  risk_fat_orangasli$AMEF |> mutate(group = "Orang Asli"),
  risk_fat_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_fat_AMEF

# Generate paired AERF + AMEF panels
risk_fat_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_fat_AERF, 
  dat_AMEF   = risk_fat_AMEF,
  outcome    = "Body Fat",
  units_AERF = "%",
  units_AMEF = "%",
  legend_pos = c(0.25, 0.85)
 )

# save
# ggsave(plot=risk_fat_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_fat_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_fat_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_fat_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 4. CRP vs INDUSTRIALIZATION
# ==============================================================================

# Combine per-group and difference AERF data (CRP on log scale)
bind_rows(
  risk_crp_turkana$AERF |> mutate(group = "Turkana"),
  risk_crp_orangasli$AERF |> 
    mutate(group = "Orang Asli"),
  risk_crp_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_crp_AERF

# Combine AMEF data. Per-group slopes are on the log scale, so back-transform
# with (exp(x) - 1) * 100 to get percent change per unit of industrialization.
# The difference AMEF is already on the percent-change scale (from
# create_pct_diff()), so it only needs renaming and class stripping.
bind_rows(
  risk_crp_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ (exp(.x) -1) * 100)),
  risk_crp_orangasli$AMEF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ (exp(.x) -1) * 100)),
  risk_crp_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_crp_AMEF

# Generate paired AERF + AMEF panels; pct_prefix=TRUE labels AMEF y-axis as %
risk_crp_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_crp_AERF, 
  dat_AMEF   = risk_crp_AMEF,
  outcome    = "CRP",
  units_AERF = "ln mg/L",
  units_AMEF = NULL,
  pct_prefix = TRUE,
  ylab_AERF_override = bquote(log~CRP~"(" * ln~mg/L * ")"),
  legend_pos = c(0.75, 0.85)
 )

# save
# ggsave(plot=risk_crp_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_crp_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_crp_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_crp_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 5. DAILY STEP COUNTS vs INDUSTRIALIZATION
# ==============================================================================

# Combine per-group and difference AERF data
bind_rows(
  risk_steps_turkana$AERF |> mutate(group = "Turkana"),
  risk_steps_orangasli$AERF |> mutate(group = "Orang Asli"),
  risk_steps_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_steps_AERF

# Combine per-group and difference AMEF data
bind_rows(
  risk_steps_turkana$AMEF |> mutate(group = "Turkana"),
  risk_steps_orangasli$AMEF |> mutate(group = "Orang Asli"),
  risk_steps_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_steps_AMEF

# Generate paired AERF + AMEF panels
risk_steps_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_steps_AERF, 
  dat_AMEF   = risk_steps_AMEF,
  outcome    = "Mean Daily Steps",
  units_AERF = NULL,
  units_AMEF = NULL,
  legend_pos = c(0.75, 0.85)
  )

# save
# ggsave(plot=risk_steps_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_steps_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_steps_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_steps_plot_diff.pdf"),
#        width = 2.5, height = 2.5)


# ==============================================================================
# 6. MVPA vs INDUSTRIALIZATION
# ==============================================================================

# Combine per-group and difference AERF data
bind_rows(
  risk_mvpa_turkana$AERF |> mutate(group = "Turkana"),
  risk_mvpa_orangasli$AERF |> mutate(group = "Orang Asli"),
  risk_mvpa_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_mvpa_AERF

# Combine per-group and difference AMEF data
bind_rows(
  risk_mvpa_turkana$AMEF |> mutate(group = "Turkana"),
  risk_mvpa_orangasli$AMEF |> mutate(group = "Orang Asli"),
  risk_mvpa_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
risk_mvpa_AMEF

# Generate paired AERF + AMEF panels
risk_mvpa_plots <- plot_AERF_AMEF(
  dat_AERF   = risk_mvpa_AERF, 
  dat_AMEF   = risk_mvpa_AMEF,
  outcome    = "MVPA",
  units_AERF = "min/day",
  units_AMEF = "min/day",
  legend_pos = c(0.75, 0.85)
 )

# save
# ggsave(plot=risk_mvpa_plots$group_plot,
#        file=file.path(base_dir, "figures/risk_mvpa_plot_groups.pdf"),
#        width = 2.5, height = 2.5)
# 
# ggsave(plot=risk_mvpa_plots$diff_plot,
#        file=file.path(base_dir, "figures/risk_mvpa_plot_diff.pdf"),
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
p1 <- tighten_pair(risk_bmi_plots$group_plot)
p2 <- tighten_pair(risk_waist_plots$group_plot)
p3 <- tighten_pair(risk_fat_plots$group_plot)
p4 <- tighten_pair(risk_crp_plots$group_plot)
p5 <- tighten_pair(risk_steps_plots$group_plot)
p6 <- tighten_pair(risk_mvpa_plots$group_plot)

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
       file=file.path(base_dir, "figures/risk_plots_groups.pdf"),
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
p1 <- tighten_pair(risk_bmi_plots$diff_plot)
p2 <- tighten_pair(risk_waist_plots$diff_plot)
p3 <- tighten_pair(risk_fat_plots$diff_plot)
p4 <- tighten_pair(risk_crp_plots$diff_plot)
p5 <- tighten_pair(risk_steps_plots$diff_plot)
p6 <- tighten_pair(risk_mvpa_plots$diff_plot)

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
       file=file.path(base_dir, "figures/risk_plots_diffs.pdf"),
       width = 5, height = 8, units = "in")

