# Analysis Code: Low Back Pain and Industrialization

Analysis code for the study of musculoskeletal pain prevalence and risk factors across four populations (Turkana, Orang Asli, Tsimane, Raramuri/Tarahumara).

## Data availability

The datasets required to run these scripts are archived separately with restricted access. Place the following CSV files in the `data/` directory:

- `Turkana_pain.csv`
- `OA_HeLP_pain.csv`
- `Tsimane_back_pain.csv`
- `Tarahumara_low_back_pain.csv`

Data are archived separately with restricted access. See the associated publication or the Zenodo metadata for this repository for the data DOI.

## Setup

1. Open `LBP_analysis.Rproj` in RStudio (or set your working directory to this folder)
2. Install the `pacman` package if not already installed: `install.packages("pacman")`
3. All other packages will be installed automatically via `pacman::p_load()` when scripts are sourced

## Scripts

Scripts are numbered in order of execution:

| Script | Description |
|--------|-------------|
| `0a_functions.R` | Shared function library (themes, descriptive helpers, causal inference, prevalence estimation) |
| `0b_setup.R` | Environment setup: loads packages, sources functions, loads and prepares datasets |
| `01_descriptive_summary.R` | Basic descriptive statistics across all four populations |
| `02_descriptive_orangasli_turkana.R` | GAM-based age-adjusted prevalence for Turkana and Orang Asli |
| `03_descriptive_tsimane.R` | GAM-based age-adjusted prevalence for Tsimane |
| `04_descriptive_tarahumara.R` | GAM-based age-adjusted prevalence for Tarahumara |
| `05_pain_industrialization.R` | Causal effect of industrialization on LBP (Turkana, Orang Asli, difference) |
| `06_risk_industrialization_turkana.R` | Causal effect of industrialization on risk factors (Turkana) |
| `07_risk_industrialization_orang_asli.R` | Causal effect of industrialization on risk factors (Orang Asli) |
| `08_risk_industrialization_diff.R` | Between-group difference in industrialization effects on risk factors |
| `09_pain_risk_turkana.R` | Causal effect of risk factors on LBP (Turkana) |
| `10_pain_risk_orang_asli.R` | Causal effect of risk factors on LBP (Orang Asli) |
| `11_pain_risk_diff.R` | Between-group difference in risk factor effects on LBP |
| `12_risk_vs_industrialization.plots.R` | Visualization: risk factors vs industrialization |
| `13_pain_vs_risk_plots.R` | Visualization: LBP vs risk factors |
| `14_pain_vs_industrialization_plots.R` | Visualization: LBP vs industrialization |

Each script sources `0b_setup.R` at the top, which handles all package loading and data preparation. Output figures are saved to `figures/` and model objects to `models/` (both created automatically).

## Statistical methods

- **Descriptive analyses**: Generalized additive models (GAMs) with sex-specific penalized thin-plate splines for age-adjusted prevalence estimation and age-prevalence curves
- **Causal analyses**: Continuous exposure effect estimation using entropy/optimal balancing weights (WeightIt), natural spline outcome models, G-computation (marginaleffects), and clarify-based simulation for between-group differences

## License

MIT
