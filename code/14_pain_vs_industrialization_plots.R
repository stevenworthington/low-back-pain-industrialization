################################################################################
#
# Pr(LBP) vs INDUSTRIALIZATION INDEX -- VISUALIZATION
# Low back pain causal analysis pipeline
#
# This script visualizes the causal effect of the industrialization index on
# the probability of low back pain (LBP). Unlike scripts 9 and 10, there is
# only a single outcome here (Pr(LBP)), so no multi-panel grid assembly is
# needed.
#
# The script:
#   - Loads all saved model .RData files from models/
#   - Combines per-group AERF data (Turkana, Orang Asli; x100 for %) with
#     the between-group difference AERF into a single tibble
#   - Repeats the same for AMEF data
#   - Calls plot_AERF_AMEF() with outcome = quote(Pr(LBP)) and
#     exposure = "industrialization_index"
#   - Saves the per-group panel to figures/pain_vs_index_plot_groups.pdf
#   - Saves the difference panel to figures/pain_vs_index_plot_diff.pdf
#
# Dependencies:
#   00b_setup.R (packages, helpers, data); saved model objects in models/
#   Packages: tidyverse, patchwork, vctrs, ggplot2
#   Helper functions: plot_AERF_AMEF()
#
# Outputs:
#   figures/pain_vs_index_plot_groups.pdf  -- 2.5 x 2.5 in (per-group panel)
#   figures/pain_vs_index_plot_diff.pdf    -- 2.5 x 2.5 in (difference panel)
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
# 1. Pr(LBP) vs INDUSTRIALIZATION INDEX
# ==============================================================================

# Combine per-group and difference AERF data. Per-group estimates are on the
# probability scale, so multiply by 100 to express as percentages. Difference
# data from clarify uses different column names, requiring rename + vec_data().
bind_rows(
  pain_lowback_turkana$AERF |> 
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_lowback_orangasli$AERF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_lowback_diff$AERF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_lowback_AERF

# Combine per-group and difference AMEF data (x100 for percent)
bind_rows(
  pain_lowback_turkana$AMEF |>
    mutate(group = "Turkana",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_lowback_orangasli$AMEF |> 
    mutate(group = "Orang Asli",
           across(c(estimate, conf.low, conf.high), ~ .x * 100)),
  pain_lowback_diff$AMEF_diff |> 
    rename(estimate = Estimate,
           conf.low = `2.5 %`,
           conf.high = `97.5 %`) |> 
    mutate(group = "Difference",
           across(c(estimate, conf.low, conf.high), vctrs::vec_data))
 ) ->
pain_lowback_AMEF

# Generate paired AERF + AMEF panels for Pr(LBP) vs industrialization
pain_lowback_plots <- plot_AERF_AMEF(
  dat_AERF   = pain_lowback_AERF, 
  dat_AMEF   = pain_lowback_AMEF,
  outcome    = quote(Pr(LBP)),
  units_AERF = NULL,
  units_AMEF = "pp/unit",
  legend_pos = c(0.25, 0.9),
  exposure   = "industrialization_index",
  xlab       = "Industrialization Index"
)

# Save per-group and difference panels as individual PDFs
ggsave(plot=pain_lowback_plots$group_plot,
       file=file.path(base_dir, "figures/pain_vs_index_plot_groups.pdf"),
       width = 2.5, height = 2.5)

ggsave(plot=pain_lowback_plots$diff_plot,
       file=file.path(base_dir, "figures/pain_vs_index_plot_diff.pdf"),
       width = 2.5, height = 2.5)


