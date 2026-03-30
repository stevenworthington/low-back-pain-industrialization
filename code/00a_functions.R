################################################################################
#
# SHARED FUNCTION LIBRARY
# Low back pain causal analysis pipeline
#
# This file defines all shared helper and utility functions used across the
# LBP analysis scripts. It is sourced by 00b_setup.R at the start of every
# analysis session.
#
# Functions are organized into the following categories:
#
#   1. Custom ggplot themes
#      - theme_bw_lbp()          : black-and-white theme for manuscript figures
#
#   2. Descriptive / data-wrangling helpers
#      - summarise_binary_outcomes() : raw prevalence (%) for binary outcomes
#      - nice_labels()               : clean variable names into readable labels
#      - count_complete_cases()      : complete-case counts per outcome + covariates
#      - summarise_numeric()         : five-number summary + NA count by group
#      - format_thousands()          : numeric axis labels with "k" suffix
#      - make_y_lab()               : plotmath y-axis labels with units
#
#   3. Multi-group descriptive prevalence estimation
#      - common_age_support()               : overlapping age range across groups
#      - fix_pr_p()                         : correct PR p-values for H0: PR = 1
#      - fit_predict_one()                  : multi-group prevalence workhorse
#                                             (GAM, age-standardized, sex PRs/diffs,
#                                              between-group PRs via delta method)
#      - fit_predict_all()                  : wrapper across multiple outcomes
#      - derivative_plateau_slopes()        : derivative/plateau analysis
#      - test_group_x_age()                 : group x age interaction test (LRT)
#      - prevalence_by_age_x_sex_curve()    : age-prevalence curve data builder
#
#   4. Causal inference (continuous exposure effects)
#      - cont_treatment_effect()     : single-group AERF/AMEF via entropy/
#                                      optimal balancing weights (WeightIt),
#                                      natural splines, G-computation
#                                      (avg_predictions), and avg_slopes
#                                      (marginaleffects)
#      - plot_AERF_AMEF()            : paired AERF + AMEF plots for per-group
#                                      and difference panels
#      - create_pct_diff()           : AMEF difference (Turkana - Orang Asli)
#                                      for logged outcomes via clarify objects
#      - cont_treatment_effect_diff(): between-group causal analysis using
#                                      clarify simulation (sim_adrf / sim)
#
# Package dependencies:
#   ggplot2, scales, patchwork, ggokabeito, cowplot          (plotting)
#   mgcv                                                     (GAMs)
#   marginaleffects                                          (avg_predictions,
#                                                             avg_comparisons,
#                                                             avg_slopes)
#   WeightIt, cobalt                                         (balancing weights,
#                                                             balance diagnostics)
#   clarify                                                  (simulation-based
#                                                             inference for diffs)
#   sandwich                                                 (cluster-robust SEs)
#   splines                                                  (natural splines in
#                                                             outcome models)
#   dplyr, tidyr, purrr, stringr, tibble                     (data wrangling)
#
################################################################################


# ==============================================================================
# 1. CUSTOM GGPLOT THEMES
# ==============================================================================

theme_bw_lbp <- function(){ 
  theme_bw() %+replace%    
    theme(plot.title        = element_text(size = 6, hjust = 0.5, 
                                    margin = margin(b = 2)),
          axis.title        = element_text(size = 7),
          axis.text         = element_text(size = 6),
          axis.title.y      = element_text(angle  = 90,          
                                           vjust  = 0.5,         
                                           margin = margin(r = 1)),
          axis.ticks        = element_line(linewidth = 0.3),
          panel.grid.minor  = element_blank(),
          panel.grid.major  = element_line(linewidth=0.2),
          strip.text        = element_text(size = 7, 
                                           margin = margin(b = 2), 
                                           lineheight = 1.1),      
          strip.background  = element_blank(),
          strip.placement   = "outside",
          strip.clip        = "off",   
          legend.position   = "inside",
          legend.key.size   = unit(0.6, "lines"),
          legend.text       = element_text(size = 6, 
                                           margin = margin(l = 0.5)),
          legend.title      = element_text(size = 7),
          legend.spacing.x  = unit(2, "pt"),
          legend.background = element_blank()
    )
}


################################################################################
### function to summarize pain prevalence
################################################################################

summarise_binary_outcomes <- function(
    data,
    pattern,
    exclude_pattern = NULL,
    label_fn = NULL
) {
  
  vars <- names(data) |>
    str_subset(pattern)
  
  if (!is.null(exclude_pattern)) {
    vars <- vars |>
      str_subset(exclude_pattern, negate = TRUE)
  }
  
  if (length(vars) == 0) {
    stop("no variables matched `pattern` (after exclusions).")
  }
  
  data |>
    select(all_of(vars)) |>
    pivot_longer(
      cols = everything(),
      names_to = "outcome",
      values_to = "value"
    ) |>
    group_by(outcome) |>
    summarise(
      percent = mean(value, na.rm = TRUE) * 100,
      .groups = "drop"
    ) |>
    mutate(
      percent = round(percent, 1),
      outcome = if (!is.null(label_fn)) label_fn(outcome) else outcome
    ) |>
    arrange(desc(percent)) |>
    select(outcome, percent)
}

nice_labels <- function(x) {
  x |>
    str_replace_all("_", " ") |>
    str_replace_all("\\s+", " ") |>
    str_trim() |>
    str_replace("\\bny\\b", "") |>
    str_squish()
}


################################################################################
### function to count complete cases
################################################################################

count_complete_cases <- function(df) {
  # identify columns that contain "_ny_01" plus the covariates
  target_cols <- grep("_ny_01", names(df), value = TRUE)
  covariates <- c("sex", "age_years", "industrialization_index")
  
  # check that covariates exist in the data
  missing_covs <- setdiff(covariates, names(df))
  if (length(missing_covs) > 0) {
    stop("Missing covariate columns: ", paste(missing_covs, collapse = ", "))
  }
  
  # compute number of complete cases for each target column + covariates
  result <- sapply(target_cols, function(var) {
    sum(complete.cases(df[, c(var, covariates)]))
  })
  
  # return result as a tibble
  tibble::tibble(
    variable = names(result),
    n_complete_cases = as.integer(result)
  )
}


################################################################################
### five-number summary + NA count for numeric variables by group
################################################################################

summarise_numeric <- function(..., vars) {
  # Accepts one or more named data frames (name = group label) and a character
  # vector of column names. Returns a tibble with min, Q1, median, Q3, max,
  # and n_na for each variable × group combination.
  dfs <- list(...)
  purrr::map_dfr(vars, function(v) {
    purrr::map_dfr(names(dfs), function(grp) {
      d <- dfs[[grp]]
      x <- d[[v]]
      tibble::tibble(
        group    = grp,
        variable = v,
        min      = format(round(min(x, na.rm = TRUE), 1), scientific = FALSE),
        q1       = format(round(quantile(x, 0.25, na.rm = TRUE), 1), scientific = FALSE),
        median   = format(round(median(x, na.rm = TRUE), 1), scientific = FALSE),
        q3       = format(round(quantile(x, 0.75, na.rm = TRUE), 1), scientific = FALSE),
        max      = format(round(max(x, na.rm = TRUE), 1), scientific = FALSE),
        n_na     = sum(is.na(x))
      )
    })
  })
}


################################################################################
### function for formatting axis labels in thousands
################################################################################

format_thousands <- function(x) {
  if (!is.numeric(x)) stop("Input must be numeric")
  sapply(x, function(y) {
    if (is.na(y)) {
      return(NA)
    } else if (abs(y) < 1000) {
      return(as.character(y))
    } else {
      rounded_value <- round(y / 1000, 1)
      return(paste0(rounded_value, "k"))
    }
  })
}


################################################################################
### functions for outcome labels
################################################################################

# helper to make y-axis labels
make_y_lab <- function(response, units = NULL, type = c("level", "delta")) {
  
  type <- match.arg(type)
  resp <- if (is.character(response)) as.name(response) else response
  
  # core: e.g., BMI  or  Δ BMI
  core <- if (type == "level") bquote(.(resp)) else bquote(Delta~.(resp))
  
  # no units? just return the core
  if (is.null(units) || (is.character(units) && (length(units) == 0 || units == ""))) {
    return(core)
  }
  
  # coerce units to something plotmath can render
  U <- units
  if (is.expression(U)) {
    Uexpr <- U[[1]]
  } else if (is.language(U)) {
    # already a plotmath call, e.g., quote(kg/m^2) or quote(mg/L)
    Uexpr <- U
  } else if (is.character(U)) {
    # treat character units as literal text (do NOT parse)
    U <- gsub("\\\\", "/", U)            # normalize accidental backslashes
    Uexpr <- if (U == "%") quote("%") else U  # "%" needs plotmath's percent token
  } else {
    # fallback: use as-is
    Uexpr <- U
  }
  
  bquote(.(core)~"("*.(Uexpr)*")")
}

# NEW FUNCTION
make_y_lab <- function(response, units = NULL, type = c("level", "delta"),
                       percent_prefix = FALSE,
                       response_display = NULL) {
  type <- match.arg(type)
  resp_token <- if (!is.null(response_display)) response_display else {
    if (is.character(response)) as.name(response) else response
  }
  
  core <- if (type == "level") bquote(.(resp_token)) else bquote(Delta~.(resp_token))
  
  if (percent_prefix && type == "delta") return(bquote("%"~.(core)))
  
  if (is.null(units) || (is.character(units) && (length(units) == 0 || units == ""))) return(core)
  
  U <- units
  if (is.expression(U)) {
    Uexpr <- U[[1]]
  } else if (is.language(U)) {
    Uexpr <- U
  } else if (is.character(U)) {
    U <- gsub("\\\\", "/", U)
    Uexpr <- if (U == "%") quote("%") else U
  } else {
    Uexpr <- U
  }
  bquote(.(core)~"("*.(Uexpr)*")")
}


################################################################################
### function for AERF and AMEF plots
################################################################################

plot_AERF_AMEF <- function(
    dat_AERF, 
    dat_AMEF, 
    df=dat,
    outcome,
    units_AERF,
    units_AMEF,
    pct_prefix         = FALSE,
    ylab_AERF_override = NULL,
    transform_expo     = "none",
    legend_pos         = c(0.25, 0.85),
    exposure           = "industrialization_score",
    xlab               = "Industrialization Index", 
    xlim               = c(3, 41.7) # extremes of ranges
    ) {
  
  # treat anything like Pr(...) as a probability outcome
  is_prob <- is.call(outcome) && identical(as.character(outcome[[1]]), "Pr")
  ylab_fmt <- if (is_prob) scales::label_percent(scale = 1) else format_thousands
  
  # xlab formatting for logged exposures
  xlab_fmt <- if (transform_expo == "log") {
    function(x) format_thousands(exp(x))
  } else {
    format_thousands
  }
  
  # breaks for logged CRP
  xlab_log_breaks <- if (transform_expo == "log") {
    log(c(0.1, 0.5, 1, 2, 5, 10, 20, 40))
  } else {
    waiver()
  }
  
  # limits
  xlab_lmts <- if (transform_expo == "log") log(xlim) else xlim 
  
  # create rug
  add_rug <- geom_rug(data = df,
    aes(x = if (transform_expo == "log") log(.data[[exposure]]) else .data[[exposure]], y = 0, color = group), 
    sides = "b", linewidth = 0.1, alpha = 0.7, length = grid::unit(0.04, "npc"), inherit.aes = FALSE,
    position = position_jitter(width = if (transform_expo == "log") 0.05 else 0.15, height = 0, seed = 123), na.rm = TRUE)
  
  dat_AERF |> 
    filter(group %in% c("Turkana", "Orang Asli")) |> 
  ggplot(aes(x = !!sym(exposure), color = group, fill = group)) +
    geom_line(aes(y = estimate), linewidth=0.3) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3, color = NA) +
    scale_color_okabe_ito(name="") +
    scale_fill_okabe_ito(name="") +
    coord_cartesian(xlim = xlab_lmts) +
    scale_x_continuous(expand = c(0, 0),  
                       breaks = xlab_log_breaks, labels = xlab_fmt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                       labels = ylab_fmt) +
    labs(x = xlab, 
         y = if (!is.null(ylab_AERF_override)) ylab_AERF_override
         else make_y_lab(response = outcome, units = units_AERF, type = "level"),
         title = "Average Exposure-Response Function") +
    theme_bw_lbp() +
    theme(legend.position.inside = legend_pos) ->
    p_plot
  
  dat_AMEF |>   
    filter(group %in% c("Turkana", "Orang Asli")) |> 
  ggplot(aes(x = !!sym(exposure), color = group, fill = group)) +
    geom_line(aes(y = estimate), linewidth=0.3) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3, color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "grey30", linewidth = 0.3) +
    add_rug +
    scale_color_okabe_ito(name="") +
    scale_fill_okabe_ito(name="") +
    coord_cartesian(xlim = xlab_lmts,
                    ylim = if (pct_prefix) c(-50, 150) else NULL) +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = xlab_log_breaks, labels = xlab_fmt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                       labels = format_thousands) +
    labs(x = xlab, 
         y = make_y_lab(response = outcome, units = units_AMEF, type="delta", 
                        percent_prefix = pct_prefix),
         title = "Average Marginal Effects Function") +
    theme_bw_lbp() +
    theme(legend.position = "none") ->
    s_plot
  
  # Combine AERF and AMEF plots
  axis_adjust <- theme(axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank())
  combined_plot <- (p_plot + axis_adjust) / 
                    plot_spacer() / 
                    s_plot + 
                    plot_layout(heights = c(1, -0.27, 1))
  
  dat_AERF |> 
    filter(group == "Difference") |> 
  ggplot(aes(x = !!sym(exposure))) +
    geom_line(aes(y = estimate), linewidth=0.3, color = "grey20") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                alpha = .3, color = NA, fill = "grey40") +
    coord_cartesian(xlim = xlab_lmts) +
    scale_x_continuous(expand = c(0, 0),  
                       breaks = xlab_log_breaks, labels = xlab_fmt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                       labels = ylab_fmt) +
    labs(x = xlab, 
         y = if (!is.null(ylab_AERF_override)) ylab_AERF_override
         else make_y_lab(response = outcome, units = units_AERF, type = "level"),
         title = "Average Exposure Response Function\n(Turkana - Orang Asli)") +
    theme_bw_lbp() -> 
    p_plot_diff
  
  dat_AMEF |>   
    filter(group == "Difference") |>  
  ggplot(aes(x = !!sym(exposure))) +
    geom_line(aes(y = estimate), linewidth=0.3, color = "grey20") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                alpha = .3, color = NA, fill = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "grey30", linewidth = 0.3) +
    coord_cartesian(xlim = xlab_lmts,
                    ylim = if (pct_prefix) c(-50, 150) else NULL) +
    scale_x_continuous(expand = c(0, 0),  
                       breaks = xlab_log_breaks, labels = xlab_fmt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                       labels = format_thousands) +
    labs(x = xlab, 
         y = make_y_lab(response = outcome, units = units_AMEF, type="delta",
                        percent_prefix = pct_prefix),
         title = "Average Marginal Effects Function\n(Turkana - Orang Asli)") +
    theme_bw_lbp() ->
    s_plot_diff
  
  # Combine AERF and AMEF plots
  axis_adjust <- theme(axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank())
  combined_plot_diff <- (p_plot_diff + axis_adjust) / 
                         plot_spacer() / 
                         s_plot_diff + 
                         plot_layout(heights = c(1, -0.3 ,1))
  
  list(group_plot=combined_plot, diff_plot=combined_plot_diff)
}


################################################################################
### function for continuous treatment causal analysis
################################################################################
cont_treatment_effect <- function(dat, treatment, outcome, covariates,
                                  outcome_type="continuous", transform_out="none", 
                                  transform_expo="none",
                                  n_moments=2, wint = TRUE, interactions = TRUE,
                                  w_meth = "energy",
                                  LQ=0.01, UQ=0.99, 
                                  xlim = c(3, 41.7), # extremes of ranges
                                  xlab = "Industrialization") {
  
  # Get treatment name
  treatment_name <- gsub("_", " ", treatment) |> tools::toTitleCase()
  
  # Ensure the treatment, outcome, and covariates are in the dat
  dat <- dat |> 
    drop_na(!!sym(treatment), !!sym(outcome), all_of(covariates))
  
  # log outcome?
  if (transform_out == "log") {
    dat[[outcome]] <- log(dat[[outcome]])
  }
  
  # log exposure?
  if (transform_expo == "log") {
    dat[[treatment]] <- log(dat[[treatment]])
  }
  
  # Fit the treatment model
  formula_treat <- formula(paste(treatment, "~", paste(covariates, collapse = " + ")))
  W <- weightit(formula = formula_treat, data = dat, moments = n_moments, int = wint, method = w_meth)
  
  # Assess balance on covariates
  balance_summary <- summary(W)
  bal_table <- bal.tab(W, poly=2)
  bal_plots <- lapply(covariates, function(var) bal.plot(W, var.name = var))
  love_plot <- love.plot(W, var.order = "unadjusted",  abs = TRUE)
  
  # Fit the outcome model
  if (interactions == TRUE) {
      formula_outcome <- formula(paste(outcome, "~ splines::ns(", treatment, ", df = 4) * (", paste(covariates, collapse = " + "), ")"))
  } else {
    formula_outcome <- formula(paste(outcome, "~ splines::ns(", treatment, ", df = 4) + (", paste(covariates, collapse = " + "), ")"))
  }
  
  if (outcome_type == "continuous") {
    fit <- lm_weightit(formula = formula_outcome, 
                        data = dat, weightit = W)
  } else if(outcome_type == "binary") {
    fit <- glm_weightit(formula = formula_outcome, 
                        family = binomial(link="logit"), 
                        data = dat, br = TRUE, weightit = W,
                        control = list(maxit=10000, type="MPL_Jeffreys"))
  } else {
    stop("Unsupported outcome type. Must be 'continuous' or 'binary'.")
  }
  
  # Representative values of treatment
  values <- seq(quantile(dat[[treatment]], probs = LQ, na.rm = TRUE), 
                quantile(dat[[treatment]], probs = UQ, na.rm = TRUE), 
                length.out = 101)
  
  # G-computation (predictions at representative values)
  variables <- list()
  variables[[treatment]] <- values
  
  if (outcome_type == "continuous") {
    p <- avg_predictions(model = fit, variables = variables, type = "response")
  } else if(outcome_type == "binary") {
    p <- avg_predictions(model = fit, variables = variables, type = "link")
    z <- qnorm(0.975)
    p <- p |> dplyr::mutate(
      eta = estimate,         # save link-scale mean
      estimate = plogis(eta), # back-transform the mean
      conf.low  = plogis(eta - z * std.error),
      conf.high = plogis(eta + z * std.error)) |>
      dplyr::select(-eta)
  } else {
    stop("Unsupported outcome type. Must be 'continuous' or 'binary'.")
  }
  
  # Plot the average exposure-response function (AERF) with tertile background
  p_plot <- ggplot(p, aes(x = !!sym(treatment))) +
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    scale_x_continuous(expand = c(0, 0), limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(x = xlab, y = "E[Y|A]") +
    theme_bw() +
    theme(legend.position = "none")
  if (transform_out == "log") {
     p_plot <- p_plot + labs(x = xlab, y = "ln E[Y|A]")
  }
  if (transform_expo == "log") {
    p_plot <- p_plot + labs(x = paste("ln", xlab), y = "E[Y|A]")
  }
  
  # Estimate the point-wise derivatives at representative values of treatment
  args <- list(grid_type = "counterfactual", newdata = dat)
  args[[treatment]] <- values
  # Use do.call to pass named arguments dynamically to datgrid
  newdat_s <- do.call(datagrid, args) # 101 * 778 = 78578
  s <- avg_slopes(model = fit, variables = treatment, newdata = newdat_s, 
                  by = treatment, slope = "dydx", type = "response")
 
  
  # Plot the average marginal effect function (AMEF)
  s_plot <- ggplot(s, aes(x = !!sym(treatment))) +
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_rug(data=dat, aes(x = !!sym(treatment), y = 0), 
             sides="b", color="grey40", inherit.aes = FALSE,
             position = position_jitter(width = 0.15, height = 0, 
                                        seed = 123)) +
    scale_x_continuous(expand = c(0, 0), limits = xlim) +
    labs(x = xlab, y = "dE[Y|A] / dA") +
    theme_bw() +
    theme(legend.position = "none")
  if (transform_out == "log") {
    s_plot <- s_plot + labs(x = xlab, y = "d ln E[Y|A] / dA")
  }
  if (transform_expo == "log") {
    s_plot <- s_plot + labs(x = paste("ln", xlab), y = "dE[Y|A] / dA")
  }
  
  # Combine AERF and AMEF plots
  axis_adjust <- theme(axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank())
  combined_plot <- (p_plot + axis_adjust) / plot_spacer() / s_plot + plot_layout(heights = c(1, -0.1 ,1))
  
  # Calculate overall mean of AMEF
  overall_amef <- avg_slopes(model = fit, variables = treatment, newdata = newdat_s, 
                             by = TRUE, slope = "dydx", type = "response")
  
  # Return a list of results
  list(
    balance_summary = balance_summary,
    bal_table = bal_table,
    bal_plots = bal_plots,
    love_plot = love_plot,
    combined_plot = combined_plot,
    overall_amef = overall_amef,
    AERF = p,
    AMEF = s
  )
}

 
################################################################################
### helper function for continuous treatment causal analysis SIMULATION
################################################################################
 
# Create AMEF "Difference" object (Turkana − Orang Asli).
# `transform_out` should be the SAME string your parent function uses: "log" or "none".
create_pct_diff <- function(amef1, amef2, transform_out) {
  transform_out <- match.arg(transform_out, c("none", "log"))
  
  # checks
  at1 <- attr(amef1, "at"); at2 <- attr(amef2, "at")
  if (!identical(at1, at2)) stop("'amef1' and 'amef2' must share the same 'at' grid.")
  var1 <- attr(amef1, "var"); var2 <- attr(amef2, "var")
  if (!identical(var1, var2)) stop("'amef1' and 'amef2' must share the same 'var'.")
  
  # pulls: draws x grid
  a1 <- unclass(amef1)
  a2 <- unclass(amef2)
  
  # point estimates (length = n_grid)
  orig1 <- as.numeric(attr(amef1, "original"))
  orig2 <- as.numeric(attr(amef2, "original"))
  
  # transform each group's slope, then subtract
  tfun <- if (transform_out == "log") function(x) 100 * (exp(x) - 1) else identity
  dd        <- tfun(a1) - tfun(a2)          # matrix: draws x grid
  orig_diff <- tfun(orig1) - tfun(orig2)    # vector: length = n_grid
  
  # column labels for the grid
  coln <- tryCatch(names(amef1), error = function(e) NULL)
  if (is.null(coln)) coln <- paste0("E[dY/d(", var1, ")|", at1, "]")
  
  # build the clarify_* object: set colnames BEFORE setting class
  out <- dd
  colnames(out) <- coln                     # <-- safe now (not yet clarify_*)
  attr(out, "original") <- stats::setNames(orig_diff, coln)
  attr(out, "at")       <- at1
  attr(out, "var")      <- var1
  attr(out, "contrast") <- "amef"
  
  # preserve the original class stack (often c("clarify_adrf","clarify_est"))
  class(out) <- class(amef1)
  
  out
}


################################################################################
### function for continuous treatment causal analysis SIMULATION
################################################################################
cont_treatment_effect_diff <- function(dat1, dat2, treatment, outcome, 
                                       covariates, outcome_type="continuous", n_moments=2,
                                       transform_out="none", 
                                       transform_expo="none", niter=1000,
                                       wint = TRUE, interactions = TRUE,
                                       val_seq = seq(5.4, 34.2, length.out=50), # range of overlap
                                       w_meth = "energy",
                                       xlim = c(3, 41.7), # extremes of ranges
                                       xlab = "Industrialization") {
  
  # Get treatment name
  treatment_name <- gsub("_", " ", treatment) |> tools::toTitleCase()
  
  # Ensure the treatment, outcome, and covariates are in the dat
  dat1 <- dat1 |> 
    drop_na(!!sym(treatment), !!sym(outcome), all_of(covariates))
  dat2 <- dat2 |> 
    drop_na(!!sym(treatment), !!sym(outcome), all_of(covariates))
  
  # log outcome?
  if (transform_out == "log") {
    dat1[[outcome]] <- log(dat1[[outcome]])
    dat2[[outcome]] <- log(dat2[[outcome]])
  }
  
  # log exposure?
  if (transform_expo == "log") {
    dat1[[treatment]] <- log(dat1[[treatment]])
    dat2[[treatment]] <- log(dat2[[treatment]])
  }
  
  # Fit the treatment models
  formula_treat <- formula(paste(treatment, "~", paste(covariates, collapse = " + ")))
  W1 <- weightit(formula = formula_treat, data = dat1, moments = n_moments, int = wint, method = w_meth)
  W2 <- weightit(formula = formula_treat, data = dat2, moments = n_moments, int = wint, method = w_meth)
  
  # Assess balance on covariates
  balance_summary1 <- summary(W1)
  bal_table1 <- bal.tab(W1, poly=2)
  bal_plots1 <- lapply(covariates, function(var) bal.plot(W1, var.name = var))
  love_plot1 <- love.plot(W1, var.order = "unadjusted",  abs = TRUE)
  balance_summary2 <- summary(W2)
  bal_table2 <- bal.tab(W2, poly=2)
  bal_plots2 <- lapply(covariates, function(var) bal.plot(W2, var.name = var))
  love_plot2 <- love.plot(W2, var.order = "unadjusted",  abs = TRUE)
  
  # Fit the outcome models
  if (interactions == TRUE) {
      formula_outcome <- formula(paste(outcome, "~ splines::ns(", treatment, ", df = 4) * (", paste(covariates, collapse = " + "), ")"))
  } else {
    formula_outcome <- formula(paste(outcome, "~ splines::ns(", treatment, ", df = 4) + (", paste(covariates, collapse = " + "), ")"))
  }
  
  if (outcome_type == "continuous") {
    fit1 <- lm_weightit(formula = formula_outcome, data = dat1, weightit = W1)
    fit2 <- lm_weightit(formula = formula_outcome, data = dat2, weightit = W2)
  } else if(outcome_type == "binary") {
    fit1 <- glm_weightit(formula = formula_outcome, family = binomial(link="logit"), 
                         data = dat1, br = TRUE, weightit = W1, control = list(maxit=10000, type="MPL_Jeffreys"))
    fit2 <- glm_weightit(formula = formula_outcome, family = binomial(link="logit"), 
                         data = dat2, br = TRUE, weightit = W2, control = list(maxit=10000, type="MPL_Jeffreys"))
  } else {
    stop("Unsupported outcome type. Must be 'continuous' or 'binary'.")
  }
  
  # Clarify
  set.seed(123)
  cl1 <- sim(fit = fit1, n = niter)
  aerf1 <- sim_adrf(sim=cl1, var=treatment, contrast="adrf", 
                    at=val_seq, n=50, type="response")
  amef1 <- sim_adrf(sim=cl1, var=treatment, contrast="amef", 
                    at=val_seq, n=50, type="response")
  cl2 <- sim(fit = fit2, n = niter)
  aerf2 <- sim_adrf(sim=cl2, var=treatment, contrast="adrf", 
                    at=val_seq, n=50, type="response")
  amef2 <- sim_adrf(sim=cl2, var=treatment, contrast="amef", 
                    at=val_seq, n=50, type="response")
  aerf_diff <- aerf1 - aerf2
  # AMEF difference: absolute % change per unit (Turkana − Orang Asli)
  if (transform_out == "log") {
    amef_diff <- create_pct_diff(amef1, amef2, transform_out = transform_out)
  } else {
    amef_diff <- amef1 - amef2  # plain difference when outcome not logged
  }
  
  aerf_diff_df <- summary(aerf_diff) |> as_tibble(rownames = treatment) |> 
    mutate(!!sym(treatment) := str_extract(!!sym(treatment), "-?[0-9.]+") |> as.numeric()) |> 
    mutate(across(where(~ inherits(., "smmry.c_")), ~ as.numeric(unclass(.))))
  amef_diff_df <- summary(amef_diff) |> as_tibble(rownames = treatment) |> 
    mutate(!!sym(treatment) := str_extract(!!sym(treatment), "-?[0-9.]+") |> as.numeric()) |> 
    mutate(across(where(~ inherits(., "smmry.c_")), ~ as.numeric(unclass(.))))
  
  # Plot the average exposure-response function (AERF) with tertile background
  p_plot <- ggplot(aerf_diff_df, aes(x = !!sym(treatment))) +
    geom_line(aes(y = Estimate)) +
    geom_ribbon(aes(ymin = `2.5 %`, ymax = `97.5 %`), alpha = .3) +
    scale_x_continuous(expand = c(0, 0), limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(x = xlab, y = "Difference in E[Y|A]") +
    theme_bw() +
    theme(legend.position = "none")
  if (transform_out == "log") {
    p_plot <- p_plot + labs(x = xlab, y = "ln Difference in E[Y|A]")
  }
  if (transform_expo == "log") {
    p_plot <- p_plot + labs(x = paste("ln", xlab), y = "Difference in E[Y|A]")
  }
  
  # Plot the average marginal effect function (AMEF)
  s_plot <- ggplot(amef_diff_df, aes(x = !!sym(treatment))) +
    geom_line(aes(y = Estimate)) +
    geom_ribbon(aes(ymin = `2.5 %`, ymax = `97.5 %`), alpha = .3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_rug(data=dat1, aes(x = !!sym(treatment), y = 0), 
             sides="b", color="grey30", inherit.aes = FALSE,
             position = position_jitter(width = 0.15, height = 0,
                                        seed = 123)) +
    geom_rug(data=dat2, aes(x = !!sym(treatment), y = 0), 
             sides="b", color="grey60", inherit.aes = FALSE,
             position = position_jitter(width = 0.15, height = 0, 
                                        seed = 123)) +
    scale_x_continuous(expand = c(0, 0), limits = xlim) +
    labs(x = xlab, y = "Difference in dE[Y|A] / dA") +
    theme_bw() +
    theme(legend.position = "none")
  if (transform_out == "log") {
    s_plot <- s_plot + labs(x = xlab, y = "Difference in ln dE[Y|A] / dA")
  }
  if (transform_expo == "log") {
    s_plot <- s_plot + labs(x = paste("ln", xlab), y = "Difference in dE[Y|A] / dA")
  }
  
  # Combine AERF and AMEF plots
  axis_adjust <- theme(axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank())
  combined_plot <- (p_plot + axis_adjust) / plot_spacer() / s_plot + plot_layout(heights = c(1, -0.1 ,1))
  
  # Return a list of results
  list(
    balance_summary1 = balance_summary1,
    bal_table1 = bal_table1,
    bal_plots1 = bal_plots1,
    love_plot1 = love_plot1,
    balance_summary2 = balance_summary2,
    bal_table2 = bal_table2,
    bal_plots2 = bal_plots2,
    love_plot2 = love_plot2,
    combined_plot = combined_plot,
    AERF_diff = aerf_diff_df,
    AMEF_diff = amef_diff_df
  )
}


################################################################################
### Common age support across groups
################################################################################

# Computes the overlapping age range across all groups in the data, with a
# data-driven upper cap, and returns three useful objects:
#
#   age_ref  – vector of *observed* ages within the overlap window
#              (used as a standardization grid for marginal averaging)
#   age_grid – regular sequence at `step`-year increments within the overlap
#              (used for plotting age-prevalence curves)
#   bounds   – named list with lo, hi, and per-group ranges
#
# The upper bound is the minimum of (a) the smallest group-specific maximum
# age and (b) the `upper_quantile` of the pooled age distribution. This
# avoids extrapolating into sparse tails while staying data-driven.
#
# Arguments:
#   data            – data frame with columns `group` and `age_years`
#   upper_quantile  – quantile of pooled ages used as a soft upper cap
#                     (default 0.99; set to 1 to use the full observed range)
#   step            – increment for the regular age_grid (default 1 year)

common_age_support <- function(data, upper_quantile = 0.99, step = 1) {

  # If no group column, treat entire dataset as one group
  if (!"group" %in% names(data)) {
    data$group <- "all"
  }

  d <- data |>
    dplyr::filter(!is.na(group), !is.na(age_years))

  # Per-group min/max
  rng <- d |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      lo = min(age_years, na.rm = TRUE),
      hi = max(age_years, na.rm = TRUE),
      .groups = "drop"
    )

  # Overlap: lower bound = max of group minimums
  lo <- max(rng$lo)

  # Upper bound: data-driven cap (pooled quantile, but no higher than the

  # smallest group maximum — you can't have common support beyond what the
  # most restricted group observed)
  hi_quantile <- as.numeric(stats::quantile(d$age_years, probs = upper_quantile,
                                            na.rm = TRUE))
  hi <- min(min(rng$hi), hi_quantile)

  # Observed ages within the overlap (for standardization)
  age_ref <- d$age_years[d$age_years >= lo & d$age_years <= hi]

  # Regular grid (for curve plotting)
  age_grid <- seq(floor(lo), floor(hi), by = step)

  list(
    age_ref  = age_ref,
    age_grid = age_grid,
    bounds   = list(lo = lo, hi = hi, per_group = rng)
  )
}


################################################################################
### Helper: correct prevalence ratio p-values to test H0: PR = 1
################################################################################

# The marginaleffects package tests H0: ratio = 0 by default for ratios.
# This helper re-computes z-statistics and p-values on the log scale so
# that the null is correctly H0: PR = 1 (i.e., log(PR) = 0).
# Expects columns: estimate + (std.error OR conf.low & conf.high).

fix_pr_p <- function(df, alpha = 0.05) {
  if ("std.error" %in% names(df)) {
    # Delta-method SE on the log scale
    se_log <- df$std.error / df$estimate
    z_new  <- log(df$estimate) / se_log
  } else if (all(c("conf.low", "conf.high") %in% names(df))) {
    # Fallback: infer SE(log) from CI width on the PR scale
    zalpha <- stats::qnorm(1 - alpha / 2)
    if (any(df$conf.low <= 0, na.rm = TRUE)) {
      stop("Cannot take log of non-positive CI bounds.")
    }
    se_log <- (log(df$conf.high) - log(df$conf.low)) / (2 * zalpha)
    z_new  <- log(df$estimate) / se_log
  } else {
    stop("Need std.error or (conf.low & conf.high) to compute p-value.")
  }
  p_new <- 2 * stats::pnorm(-abs(z_new))
  df$z <- z_new
  df[["Pr(>|z|)"]] <- p_new
  df
}


################################################################################
### Multi-group prevalence estimation with within- and between-group
### prevalence ratios and differences (marginaleffects-based)
################################################################################

# Main workhorse for computing prevalence estimates across groups. For each
# outcome × group, it fits a sex-specific GAM and produces several estimands:
#
#   1. Age-standardized prevalence by sex — predictions averaged over a common
#      age grid with equal weight per age point (option 2 in the standardization
#      notes). This is the primary estimand for between-group comparisons.
#   2. Pooled prevalence per group — collapsed across sex. With pool = "equal"
#      this is age- and sex-standardized; with pool = "sample" the
#      sex weights reflect the observed composition, which is equivalent to
#      age-standardized only.
#   3. Within-group sex prevalence ratio (Male/Female) via marginaleffects
#      delta method, with p-value corrected for H0: PR = 1.
#   4. Within-group sex prevalence difference (Male - Female, percentage points)
#      via marginaleffects delta method.
#   5. Between-group prevalence ratio (with CI and % change) via delta method
#      on independent pooled estimates from per-group models.
#
# The optional `vcov` argument supports cluster-robust or other sandwich
# estimators. Pass a formula (e.g., ~id) for cluster-robust SEs — this is
# converted internally to a function via sandwich::vcovCL, which is needed
# because the formula shorthand doesn't work with mgcv::gam objects. You can
# also pass a function or a matrix directly.
#
# Returns a list with: by_sex, pooled, age_pvals, sex_pr, sex_diff,
# group_contrast_pr.

fit_predict_one <- function(
    outcome_var,
    df          = dat,
    age_grid    = age_ref,
    pool        = c("equal", "sample"),
    group_order = NULL,
    vcov        = NULL
) {
  pool <- match.arg(pool)

  # --- data prep ---------------------------------------------------------------

  d <- df |>
    dplyr::filter(!is.na(.data[[outcome_var]]),
                  !is.na(group), !is.na(sex), !is.na(age_years)) |>
    dplyr::mutate(
      group     = factor(group),
      sex       = factor(sex, levels = c("Female", "Male")),
      age_years = as.numeric(age_years)
    )

  # --- per-group fits & outputs ------------------------------------------------

  per_group <- d |>
    dplyr::group_split(group) |>
    purrr::map(function(df_g) {
      if (nrow(df_g) == 0) return(NULL)

      # Keep consistent factor levels inside the split
      df_g <- df_g |>
        dplyr::mutate(
          group = factor(group, levels = levels(d$group)),
          sex   = factor(sex,   levels = levels(d$sex))
        )

      form <- stats::as.formula(
        paste0(outcome_var, " ~ sex + s(age_years, by = sex, k = 6, bs = 'tp')")
      )
      gm <- mgcv::gam(form, data = df_g, family = binomial(),
                      method = "REML", select = TRUE)

      # Resolve vcov for marginaleffects. Formula syntax (e.g., ~id) doesn't
      # work with mgcv::gam because sandwich lacks estfun/bread methods for
      # GAMs. Convert formulas to a function that computes cluster-robust vcov
      # from the group's data directly.
      if (inherits(vcov, "formula")) {
        cluster_col <- all.vars(vcov)
        local_df_g  <- df_g  # capture for closure
        me_vcov <- function(m) {
          sandwich::vcovCL(m, cluster = local_df_g[[cluster_col]], type = "HC0")
        }
      } else {
        me_vcov <- if (!is.null(vcov)) vcov else TRUE
      }

      sex_levels <- levels(df_g$sex)

      # Context for grids (drop group to avoid datagrid confusion)
      df_ctx <- df_g |>
        dplyr::select(-group) |>
        dplyr::mutate(sex = factor(sex, levels = sex_levels))

      # Standardization grid: age × sex
      nd_by_sex <- marginaleffects::datagrid(
        sex       = factor(sex_levels, levels = sex_levels),
        age_years = age_grid,
        newdata   = df_ctx
      )

      # 1) Age-standardized prevalence by sex (%)
      by_sex <- marginaleffects::avg_predictions(
        model   = gm,
        newdata = nd_by_sex,
        by      = "sex",
        type    = "response",
        vcov    = me_vcov
      ) |>
        dplyr::mutate(
          group     = as.character(df_g$group[1]),
          estimate  = estimate * 100,
          conf.low  = conf.low  * 100,
          conf.high = conf.high * 100
        ) |>
        dplyr::select(group, sex, estimate, conf.low, conf.high) |>
        tibble::as_tibble()

      # 2) Pooled prevalence per group (%), equal vs sample-weighted over sex
      #    Retains std.error (on % scale) for between-group PR delta method.
      sex_prop <- NULL
      if (pool == "sample") {
        sex_prop <- df_g |>
          dplyr::count(sex) |>
          dplyr::mutate(w = n / sum(n))
        nd_pool <- nd_by_sex |>
          dplyr::left_join(sex_prop, by = "sex") |>
          dplyr::mutate(w = dplyr::if_else(is.na(w), 0, w))

        pooled <- marginaleffects::avg_predictions(
          model   = gm,
          newdata = nd_pool,
          type    = "response",
          wts     = nd_pool$w,
          vcov    = me_vcov
        )
      } else {
        pooled <- marginaleffects::avg_predictions(
          model   = gm,
          newdata = nd_by_sex,
          type    = "response",
          vcov    = me_vcov
        )
      }

      pooled <- pooled |>
        dplyr::mutate(
          group     = as.character(df_g$group[1]),
          estimate  = estimate  * 100,
          std.error = std.error * 100,
          conf.low  = conf.low  * 100,
          conf.high = conf.high * 100
        ) |>
        dplyr::select(group, estimate, std.error, conf.low, conf.high) |>
        tibble::as_tibble()

      # 3) Smooth p-values for the age effect by sex (from GAM summary).
      #    Tests H0: age-prevalence curve is flat for each sex.
      st <- summary(gm)$s.table
      age_pvals <- tibble::tibble(
        group = as.character(df_g$group[1]),
        sex   = c("Female", "Male"),
        p_age = st[grep("s\\(age_years\\):sex(Female|Male)", rownames(st)), "p-value"]
      )

      #    fix_pr_p() corrects the p-value to test H0: PR = 1 (not PR = 0).
      sex_pr <- marginaleffects::avg_comparisons(
        model      = gm,
        newdata    = nd_by_sex,
        variables  = "sex",
        comparison = "ratio",
        type       = "response",
        vcov       = me_vcov
      ) |>
        tibble::as_tibble() |>
        fix_pr_p() |>
        dplyr::mutate(
          group      = as.character(df_g$group[1]),
          pct_change = 100 * (estimate - 1),
          pct.low    = 100 * (conf.low - 1),
          pct.high   = 100 * (conf.high - 1)
        ) |>
        dplyr::select(group, term, contrast, estimate, conf.low, conf.high,
                      z, p.value = `Pr(>|z|)`, pct_change, pct.low, pct.high)

      # 4) Within-group sex prevalence difference (Male - Female, %)
      sex_diff <- marginaleffects::avg_comparisons(
        model      = gm,
        newdata    = nd_by_sex,
        variables  = "sex",
        comparison = "difference",
        type       = "response",
        vcov       = me_vcov
      ) |>
        tibble::as_tibble() |>
        dplyr::mutate(
          group     = as.character(df_g$group[1]),
          estimate  = estimate  * 100,
          conf.low  = conf.low  * 100,
          conf.high = conf.high * 100
        ) |>
        dplyr::select(group, term, contrast, estimate, conf.low, conf.high,
                      statistic, p.value)

      list(
        by_sex    = by_sex,
        pooled    = pooled,
        age_pvals = age_pvals,
        sex_pr    = sex_pr,
        sex_diff  = sex_diff
      )
    })

  # Bind across groups
  by_sex_tbl    <- purrr::map(per_group, "by_sex")    |> dplyr::bind_rows()
  pooled_tbl    <- purrr::map(per_group, "pooled")    |> dplyr::bind_rows()
  age_pvals_tbl <- purrr::map(per_group, "age_pvals") |> dplyr::bind_rows()
  sex_pr_tbl    <- purrr::map(per_group, "sex_pr")    |> dplyr::bind_rows()
  sex_diff_tbl  <- purrr::map(per_group, "sex_diff")  |> dplyr::bind_rows()

  # 5) Between-group PR via delta method on independent pooled estimates.
  #    Since models are fit separately per group, the estimates are independent,
  #    so Var(log(PR)) = Var(log(p_tgt)) + Var(log(p_ref)) on the log scale.
  grp_contrast_pr <- NULL
  if (!is.null(pooled_tbl) && nrow(pooled_tbl) >= 2) {
    ordered_pooled <- pooled_tbl
    if (!is.null(group_order)) {
      keep <- intersect(group_order, ordered_pooled$group)
      ordered_pooled <- ordered_pooled |>
        dplyr::filter(group %in% keep) |>
        dplyr::mutate(group = factor(group, levels = keep)) |>
        dplyr::arrange(group)
    }

    ref_row <- ordered_pooled[1, ]
    others  <- ordered_pooled[-1, , drop = FALSE]
    zcrit   <- stats::qnorm(0.975)

    grp_contrast_pr <- purrr::map_dfr(seq_len(nrow(others)), function(i) {
      tgt <- others[i, ]
      pr_est  <- tgt$estimate / ref_row$estimate
      # SE on log scale: delta method with independence of separate models
      se_log  <- sqrt((tgt$std.error / tgt$estimate)^2 +
                      (ref_row$std.error / ref_row$estimate)^2)
      z_val   <- log(pr_est) / se_log
      p_val   <- 2 * stats::pnorm(-abs(z_val))
      ci_lo   <- exp(log(pr_est) - zcrit * se_log)
      ci_hi   <- exp(log(pr_est) + zcrit * se_log)

      tibble::tibble(
        term       = "group",
        contrast   = paste0(tgt$group, " / ", ref_row$group),
        estimate   = pr_est,
        conf.low   = ci_lo,
        conf.high  = ci_hi,
        z          = z_val,
        p.value    = p_val,
        pct_change = 100 * (pr_est - 1),
        pct.low    = 100 * (ci_lo - 1),
        pct.high   = 100 * (ci_hi - 1)
      )
    })
  }

  list(
    by_sex            = by_sex_tbl,
    pooled            = pooled_tbl,
    age_pvals         = age_pvals_tbl,
    sex_pr            = sex_pr_tbl,
    sex_diff          = sex_diff_tbl,
    group_contrast_pr = grp_contrast_pr
  )
}


################################################################################
### Wrapper: run fit_predict_one() across multiple outcomes and consolidate
################################################################################

# Calls fit_predict_one() for each outcome in `outcomes`, then stacks each
# result component into its own tidy tibble with a `location` column added
# via `lab_map`.
#
# Arguments:
#   outcomes    – character vector of binary outcome column names
#   lab_map     – named character vector mapping outcome names to pretty labels
#   ...         – additional arguments passed to fit_predict_one()
#                 (df, age_grid, pool, group_order, vcov)
#
# Returns a named list of seven tibbles, each with a `location` column:
#   by_sex, pooled, age_pvals, sex_pr, sex_diff, group_contrast_pr.
# The location column is character; apply factor() in the caller as needed.

fit_predict_all <- function(outcomes, lab_map, ...) {

  raw <- purrr::imap(setNames(outcomes, outcomes), function(outcome, nm) {
    fit_predict_one(outcome, ...)
  })

  # Helper: extract one component from each outcome's result, add location label.
  # Silently skips NULLs (e.g., group_contrast_pr when only one group exists).
  bind_component <- function(component) {
    purrr::imap_dfr(raw, function(res, nm) {
      piece <- res[[component]]
      if (!is.null(piece) && nrow(piece) > 0) {
        piece |> dplyr::mutate(location = lab_map[[nm]])
      }
    })
  }

  list(
    by_sex            = bind_component("by_sex"),
    pooled            = bind_component("pooled"),
    age_pvals         = bind_component("age_pvals"),
    sex_pr            = bind_component("sex_pr"),
    sex_diff          = bind_component("sex_diff"),
    group_contrast_pr = bind_component("group_contrast_pr")
  )
}


################################################################################
### Derivative plateau analysis and monotonicity check (marginaleffects)
################################################################################

# Estimates the age at which the derivative of the age-prevalence curve
# (on the link/logit scale) is no longer significantly positive. Uses
# marginaleffects::avg_slopes to compute pointwise Wald-based derivatives
# and CIs. Also reports monotonicity: whether the derivative is significantly
# positive across the entire age grid (i.e., prevalence is strictly increasing).
#
# Arguments:
#   df_g      – data for one group (must contain outcome, sex, age_years)
#   outcome   – name of the binary outcome column
#   sex_level – which sex to evaluate ("Female" or "Male")
#   k, bs     – GAM smooth parameters
#   age_grid  – numeric vector of ages (default: 0.25-year grid over observed range)
#
# Returns a list with:
#   sex, increasing (TRUE if derivative CI lower > 0 everywhere),
#   frac_positive (fraction of grid where lower > 0),
#   min_lower, max_upper (extremes of derivative CI bounds),
#   plateau_age (upper bound of longest contiguous positive-derivative run),
#   positive_range (age span of that run), band (full derivative CI tibble).

derivative_plateau_slopes <- function(
    df_g,
    outcome   = "low_back_pain_ny_01",
    sex_level = c("Female", "Male"),
    k = 6, bs = "tp",
    age_grid  = NULL
) {
  sex_level <- match.arg(sex_level)
  df_g <- droplevels(df_g)

  # Fit per-group GAM with sex-specific smooths
  form <- stats::as.formula(
    paste0(outcome, " ~ sex + s(age_years, by = sex, k = ", k, ", bs = '", bs, "')")
  )
  gm <- mgcv::gam(form, data = df_g, family = binomial(), method = "REML", select = TRUE)

  # Dense age grid (0.25-year increments by default)
  if (is.null(age_grid)) {
    r <- range(df_g$age_years, na.rm = TRUE)
    age_grid <- seq(r[1], r[2], by = 0.25)
  }

  # Prediction grid for one sex
  nd <- marginaleffects::datagrid(
    age_years = age_grid,
    sex       = factor(sex_level, levels = levels(df_g$sex)),
    newdata   = df_g |> dplyr::select(-group)
  )

  # Pointwise derivatives d(eta)/d(age) with Wald CIs on the logit scale
  sl <- marginaleffects::avg_slopes(
    model     = gm,
    variables = "age_years",
    newdata   = nd,
    by        = "age_years",
    type      = "link"
  ) |>
    dplyr::transmute(
      age      = .data$age_years,
      lower    = .data$conf.low,
      upper    = .data$conf.high,
      est      = .data$estimate,
      positive = lower > 0
    )

  band <- tibble::tibble(
    age      = sl$age,
    lower    = sl$lower,
    upper    = sl$upper,
    positive = sl$positive
  )

  # Find longest contiguous run where lower > 0
  rle_pos <- rle(band$positive)
  ends    <- cumsum(rle_pos$lengths)
  starts  <- ends - rle_pos$lengths + 1
  blocks  <- tibble::tibble(start = starts, end = ends, val = rle_pos$values) |>
    dplyr::filter(val) |>
    dplyr::mutate(len = end - start + 1) |>
    dplyr::arrange(dplyr::desc(len))

  # Monotonicity summary: is the derivative significantly > 0 everywhere?
  increasing  <- all(band$positive)
  frac_pos    <- mean(band$positive)
  min_lower   <- min(band$lower, na.rm = TRUE)
  max_upper   <- max(band$upper, na.rm = TRUE)

  if (nrow(blocks) == 0) {
    return(list(
      sex            = sex_level,
      increasing     = increasing,
      frac_positive  = frac_pos,
      min_lower      = min_lower,
      max_upper      = max_upper,
      plateau_age    = NA_real_,
      positive_range = NULL,
      band           = band
    ))
  }

  top       <- blocks[1, ]
  pos_range <- band$age[top$start:top$end]

  list(
    sex            = sex_level,
    increasing     = increasing,
    frac_positive  = frac_pos,
    min_lower      = min_lower,
    max_upper      = max_upper,
    plateau_age    = max(pos_range),
    positive_range = c(min(pos_range), max(pos_range)),
    band           = band
  )
}


################################################################################
### Group × age interaction test (nested GAM LRT)
################################################################################

# Tests whether the age-prevalence curve shape differs across groups within
# each sex, using a nested GAM likelihood ratio test. The null model has a
# single shared age smooth; the alternative has group-specific age smooths.
# Uses ML (not REML) and select=FALSE for a valid LRT.
#
# Arguments:
#   data    – data frame with outcome, group, sex, age_years columns
#   outcome – name of the binary outcome column
#   k, bs   – GAM smooth parameters
#
# Returns a tibble with one row per sex: sex, p_group_x_age.

test_group_x_age <- function(
    data,
    outcome = "low_back_pain_ny_01",
    k       = 6,
    bs      = "tp"
) {

  d <- data |>
    dplyr::filter(!is.na(.data[[outcome]]), !is.na(group), !is.na(sex), !is.na(age_years)) |>
    dplyr::mutate(
      group     = factor(group),
      sex       = factor(sex),
      age_years = as.numeric(age_years)
    )

  test_one_sex <- function(df_sex) {
    df_sex <- droplevels(df_sex)
    if (length(levels(df_sex$group)) < 2) return(NA_real_)
    df_sex <- dplyr::mutate(df_sex, grp = factor(group))

    m0 <- mgcv::gam(
      stats::as.formula(paste0(outcome, " ~ grp + s(age_years, k = ", k, ", bs = '", bs, "')")),
      data = df_sex, family = binomial(), method = "ML", select = FALSE
    )
    m1 <- mgcv::gam(
      stats::as.formula(paste0(outcome, " ~ grp + s(age_years, by = grp, k = ", k, ", bs = '", bs, "')")),
      data = df_sex, family = binomial(), method = "ML", select = FALSE
    )

    out <- try(stats::anova(m0, m1, test = "Chisq")$`Pr(>Chi)`[2], silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else out
  }

  tibble::tibble(
    sex           = c("Female", "Male"),
    p_group_x_age = c(
      d |> dplyr::filter(sex == "Female") |> test_one_sex(),
      d |> dplyr::filter(sex == "Male")   |> test_one_sex()
    )
  )
}


################################################################################
### Build age-prevalence curve data for multiple outcomes
################################################################################

# Constructs a tidy tibble of GAM-predicted prevalence (%) by age, sex, group,
# and outcome. Uses Wald CIs on the link scale, back-transformed to [0, 100].
# Designed for faceted plots showing prevalence curves for all pain locations.
#
# Arguments:
#   outcomes      – character vector of binary outcome column names
#   df            – data frame with outcome(s), group, sex, age_years
#   lab_map       – optional named vector mapping outcome names to pretty labels
#   k, bs         – GAM smooth parameters
#   upper_quantile – quantile of pooled ages for data-driven upper cap
#                    (default 0.99; set to 1 to use the full observed range)
#   age_step       – age grid increment (default 1 year)
#   vcov           – optional vcov specification for marginaleffects
#                    (e.g., ~id for cluster-robust SEs; NULL = model default)
#
# Returns: tibble with columns outcome, label, group, age_years, sex,
# estimate, conf.low, conf.high (all on percent scale).

prevalence_by_age_x_sex_curve <- function(
    outcomes,
    df,
    lab_map        = NULL,
    k = 6, bs = "tp",
    upper_quantile = 0.99,
    age_step       = 1,
    vcov           = NULL
) {
  stopifnot(is.character(outcomes), length(outcomes) >= 1)

  # Safe name lookup helper
  `%||%nm%` <- function(x, name, fallback) {
    if (!is.null(x) && !is.null(names(x)) && name %in% names(x)) x[[name]] else fallback
  }

  z <- stats::qnorm(0.975)

  purrr::map_dfr(outcomes, function(outcome) {
    d <- df |>
      dplyr::filter(!is.na(.data[[outcome]]),
                    !is.na(group), !is.na(sex), !is.na(age_years)) |>
      dplyr::mutate(
        group     = factor(group),
        sex       = factor(sex),
        age_years = as.numeric(age_years)
      )

    if (nrow(d) == 0) {
      return(tibble::tibble(
        outcome = character(), label = character(), group = factor(),
        age_years = numeric(), sex = factor(),
        estimate = numeric(), conf.low = numeric(), conf.high = numeric()
      ))
    }

    # Common age support for this outcome (per-outcome because missingness
    # patterns may differ across outcomes, changing group-specific ranges)
    cas <- common_age_support(d, upper_quantile = upper_quantile, step = age_step)
    lo       <- cas$bounds$lo
    hi       <- cas$bounds$hi
    age_grid <- cas$age_grid

    # Pretty label
    pretty_default <- tools::toTitleCase(gsub("_", " ", outcome, fixed = TRUE))
    label <- `%||%nm%`(lab_map, outcome, pretty_default)

    # Fit one GAM per group, predict on link scale, back-transform
    per_group <- d |>
      dplyr::group_split(group) |>
      purrr::map(function(df_g) {
        if (nrow(df_g) == 0) return(NULL)

        y <- df_g[[outcome]]
        if (length(unique(na.omit(y))) < 2) return(NULL)

        gm <- mgcv::gam(
          stats::as.formula(paste0(outcome, " ~ sex + s(age_years, by = sex, k = ", k, ", bs = '", bs, "')")),
          data = df_g, family = binomial(), method = "REML", select = TRUE
        )

        # Resolve vcov (see fit_predict_one for rationale)
        if (inherits(vcov, "formula")) {
          cluster_col <- all.vars(vcov)
          local_df_g  <- df_g
          me_vcov <- function(m) {
            sandwich::vcovCL(m, cluster = local_df_g[[cluster_col]], type = "HC0")
          }
        } else {
          me_vcov <- if (!is.null(vcov)) vcov else TRUE
        }

        marginaleffects::avg_predictions(
          model   = gm,
          newdata = marginaleffects::datagrid(
            age_years = age_grid,
            sex       = levels(df_g$sex)
          ),
          by   = c("age_years", "sex"),
          type = "link",
          vcov = me_vcov
        ) |>
          dplyr::mutate(
            .eta      = estimate,
            estimate  = plogis(.eta) * 100,
            conf.low  = plogis(.eta - z * std.error) * 100,
            conf.high = plogis(.eta + z * std.error) * 100,
            group     = df_g$group[1]
          ) |>
          dplyr::select(group, age_years, sex, estimate, conf.low, conf.high)
      }) |>
      dplyr::bind_rows()

    if (nrow(per_group) == 0) {
      return(tibble::tibble(
        outcome = character(), label = character(), group = factor(),
        age_years = numeric(), sex = factor(),
        estimate = numeric(), conf.low = numeric(), conf.high = numeric()
      ))
    }

    per_group |>
      dplyr::mutate(outcome = outcome, label = label, .before = 1)
  })
}


