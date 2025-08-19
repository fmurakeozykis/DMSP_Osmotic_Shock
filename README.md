# DMSP_Osmotic_Shock

Honours Project (Advanced Y3). Single-script analysis of osmotic shock effects on DMSP cycling in Arctic algae. Computes t0→t24 physiology deltas, fits mixed models (Type-III), runs *emmeans* contrasts (Holm), and exports publication-ready tables/figures.

---

# Physiology t0→t24 on DMS(P)t — README

This repository contains a single R script that analyzes physiology changes from **t0** to **t24** in **DMS(P)t** for two Arctic algal species. It builds paired (within-bottle) deltas, fits mixed models, computes Holm-adjusted post-hoc contrasts, and exports paper-ready tables and boxplots.

---

- **Inputs:** `clean_chr.csv`, `clean_pol_modified.csv` (or `clean_pol.csv`)
- **Outputs (in `results/`):** LaTeX + TXT **ALL contrasts**; boxplots as PNGs; one multi-page PDF
- **Model:** Δ(t24−t0) with random intercept for bottle; Type-III ANOVA; `emmeans` with Holm
- **Figures:** t24 direct metrics + ALT growth-corrected metrics

---

## What the script does (at a glance)

1. **Loads data** from two CSV files (`clean_chr.csv`, `clean_pol_modified.csv` or `clean_pol.csv`).
2. **Harmonizes factors** (`species`, `sal_adapt`→A25/A45, `sal_exp`→E25/E45, `repl`) and standardizes `time` to **t0/t24**.
3. **Pairs t0 with t24 within bottle** (`bottle = species × sal_adapt × repl`) and computes:
   - Δ mean cell size (`d_size`), Δ cell density (`d_dens`), Δ yield (Fv/Fm; `d_yield`)
   - Percent change for each Δ (`pct_size`, `pct_dens`, `pct_yield`)
4. **Fits a Δ-LMM:** delta ~ species * sal_adapt * sal_exp + (1 | bottle)

5. Falls back to `lm()` if `lmer()` fails.
5. **Runs post-hoc contrasts via `emmeans`** with **Holm** adjustment:
- Species within `sal_adapt × sal_exp`
- `E25 vs E45` within `species × sal_adapt`
- `A25 vs A45` within `species × sal_exp`
- Marginal effects and full 3-way cells (Tukey)
6. **Exports:**
- `results/ALL_contrasts_combined.tex` (LaTeX table of *all* contrasts)
- `results/ALL_contrasts_combined.txt` (human-readable summary)
- Boxplots (`results/figs/*.png`) and a multi-page PDF (`results/figs/all_boxplots_t24.pdf`)

---

## File layout & inputs

Place the script in a project folder that also contains:


### Required columns (by name)

- **Keys:** `species`, `sal_adapt` (25/45), `sal_exp` (25/45), `repl`, `time` (0/24 or t0/t24), `id`
- **For Δ-physiology:** `yield`, `cell_dens_avr`, `cell_size_avr`
- **For t24 plots / ALT metrics:**
  - From `DMSPp`: `u_poc_h_1` (μ-POC), `u_dmsp_h_1` (μ-DMSP), `u_dmsp_u_poc`, `dmsp_c_poc_mol_mol`, `total_poc_mg_l`
  - From `DMS(P)t`: `d3_p_taken_up`, `d3_p_lost_demethylated`, `dms_puptake_c_poc_mol_mol`, `dmsp_uptake_of_total_dmsp` (optional)

> The script maps `sal_adapt` → **A25/A45**, `sal_exp` → **E25/E45** and standardizes `time` to **t0/t24**.  
> Missing columns lead to **skipped plots** with an informative console message.

---

## Outputs

- **Tables**
  - `results/ALL_contrasts_combined.tex` — LaTeX table with effect sizes: `estimate [95% CI]; p (Holm)`
  - `results/ALL_contrasts_combined.txt` — same content in plain text
- **Figures**
  - `results/figs/box_<metric>.png` — one PNG per metric at t24
  - `results/figs/all_boxplots_t24.pdf` — multi-page PDF with all boxplots
- **Console**
  - QC counts of paired branches per (`species × sal_adapt × sal_exp`)
  - Messages when columns are missing or models are skipped

---

## How to run

### 1) Clone and open in R/RStudio

```bash
git clone <your-repo-url>
cd <your-repo>
install.packages(c(
  "dplyr","tidyr","readr","ggplot2","lme4","lmerTest","emmeans",
  "car","broom","xtable","rlang","tibble","purrr"
))
source("physiology_t0t24_DMSPt_ALL_contrasts.R")  # use your actual filename

Reproducibility notes

The script temporarily sets sum-to-zero contrasts and restores your original options on exit:

.old_contr <- options("contrasts")
options(contrasts = c("contr.sum","contr.poly"))
on.exit(do.call(options, .old_contr), add = TRUE)


Post-hoc P-values use Holm adjustment by default (adjust = "holm").

Helper functions safe_filename() and .bind_if() are defined if absent to avoid pipeline conflicts.

Statistical model & contrasts

Model (Δ outcome per bottle):

delta ~ species * sal_adapt * sal_exp + (1 | bottle)


where bottle = interaction(species, sal_adapt, repl)

ANOVA: car::Anova(fit, type = 3) (Type-III, sum-to-zero coding)

emmeans contrasts: species-wise, exposure-wise, acclimation-wise; plus marginal effects and Tukey on 3-way cells

Effect reporting: estimate with 95% CI and Holm-adjusted p

Figures

Direct t24 metrics: μ-POC, μ-DMSP (from DMSPp or DMS), μ-DMSP/μ-POC, DMSP-C:POC, raw D3-DMSP uptake/loss, etc.

ALT (growth-corrected) metrics:

uptake_alt = (d3_p_taken_up × μ_POC) / POC_t24

loss_alt = (d3_p_lost_demethylated × μ_POC) / POC_t24

Figures are saved under results/figs/. Optional CLD letters are implemented but off by default.

Assumptions & notes

Baseline density scaling: dens_t0_bottle = 0.6 × cell_dens_avr (documented for transparency).

Pairing: Δ is computed per bottle (unique species × sal_adapt × repl), then compared across exposure (E25 vs E45) and other factors.

Multiple testing: Holm across each contrast family; Tukey only for the full 3-way grid.

Robustness: If there are < 8 paired rows, the Δ-model for that metric is skipped with a console message.
