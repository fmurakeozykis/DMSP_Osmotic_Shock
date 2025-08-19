# ============================================================
# Physiology t0 → t24 on DMS(P)t, merged into ALL_contrasts
# ============================================================

# --- Working dir (edit if needed) ---
setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")

# --- Packages ---
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(lme4); library(lmerTest); library(emmeans); library(car)
  library(broom); library(xtable); library(rlang); library(tibble)
})

# Respect the rest of your pipeline if there are conflicts:
# Use sum-to-zero contrasts only inside this block and restore afterwards.
.old_contr <- options("contrasts")
options(contrasts = c("contr.sum","contr.poly"))
on.exit(do.call(options, .old_contr), add = TRUE)

# --- Safe helpers (defined only if missing) ---
if (!exists("safe_filename", mode = "function")) {
  safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("_+$", "", x)
    tolower(x)
  }
}
if (!exists(".bind_if", mode = "function")) {
  .bind_if <- function(lst) {
    lst <- lst[!vapply(lst, function(x) is.null(x) || !nrow(x), logical(1))]
    if (!length(lst)) return(NULL)
    dplyr::bind_rows(lst)
  }
}

# ALL-contrasts exporters (define light fallbacks if not already present)
if (!exists("make_family_table_tex_all", mode = "function")) {
  make_family_table_tex_all <- function(family_list, results, filename, sci_p_below = 1e-3){
    fmt_num <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, format="f", digits=d))
    fmt_p   <- function(p) ifelse(is.na(p), "",
                                  ifelse(p < sci_p_below, formatC(p, format="e", digits=2),
                                         formatC(p, format="f", digits=3)))
    blocks <- list()
    for (nm in family_list) {
      res <- results[[nm]]; if (is.null(res) || is.null(res$contrasts) || !nrow(res$contrasts)) next
      dfc <- res$contrasts %>%
        mutate(
          estimate = as.numeric(estimate),
          SE       = as.numeric(SE),
          lower.CL = as.numeric(lower.CL),
          upper.CL = as.numeric(upper.CL),
          p.value  = as.numeric(p.value),
          effect   = paste0(fmt_num(estimate,3), " [", fmt_num(lower.CL,3), ", ",
                            fmt_num(upper.CL,3), "]; p=", fmt_p(p.value))
        ) %>%
        select(label, contrast_set, contrast, effect)
      blocks[[length(blocks)+1]] <- dfc
    }
    out_tab <- .bind_if(blocks)
    if (is.null(out_tab)) out_tab <- tibble(label="—", contrast_set="—", contrast="No contrasts", effect="—")
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    xt <- xtable::xtable(out_tab,
                         caption = "ALL contrasts (Holm p shown; no filtering)",
                         label   = paste0("tab:", safe_filename(basename(filename))))
    out_tex <- paste0(filename, ".tex")
    print(xt, include.rownames = FALSE, file = out_tex, sanitize.text.function = identity)
    cat("Wrote LaTeX to: ", normalizePath(out_tex, winslash="/"), "\n")
    invisible(out_tab)
  }
}
if (!exists("make_family_txt_all", mode = "function")) {
  make_family_txt_all <- function(family_list, results, filename_txt, sci_p_below = 1e-3){
    fmt_num <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
    fmt_p   <- function(p) ifelse(is.na(p), "NA",
                                  ifelse(p < sci_p_below, formatC(p, format="e", digits=2),
                                         formatC(p, format="f", digits=3)))
    lines <- c("==== ALL contrasts (no filtering) ====")
    for (nm in family_list) {
      res <- results[[nm]]; if (is.null(res) || is.null(res$contrasts) || !nrow(res$contrasts)) next
      lines <- c(lines, sprintf("\n-- %s --", res$label %||% nm))
      dfc <- res$contrasts %>%
        mutate(
          estimate = as.numeric(estimate),
          lower.CL = as.numeric(lower.CL),
          upper.CL = as.numeric(upper.CL),
          p.value  = as.numeric(p.value),
          eff_str  = paste0(fmt_num(estimate,3), " [", fmt_num(lower.CL,3), ", ",
                            fmt_num(upper.CL,3), "]; p=", fmt_p(p.value))
        )
      for (cs in unique(dfc$contrast_set)) {
        lines <- c(lines, paste0("  * ", cs))
        sub <- dfc %>% filter(contrast_set == cs)
        add <- paste0("     - ", sub$contrast, ": ", sub$eff_str)
        lines <- c(lines, add)
      }
    }
    dir.create(dirname(filename_txt), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, filename_txt)
    cat("Wrote TXT to: ", normalizePath(filename_txt, winslash="/"), "\n")
    invisible(filename_txt)
  }
}

# ---------- Load & harmonise ----------
chr <- read_csv("clean_chr.csv",          show_col_types = FALSE)
pol <- read_csv("clean_pol_modified.csv", show_col_types = FALSE)

setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/statistical analysis")

df0 <- bind_rows(chr, pol) %>%
  mutate(
    time_std = case_when(
      as.character(time) %in% c("0", "t0")   ~ "t0",
      as.character(time) %in% c("24","t24")  ~ "t24",
      TRUE ~ NA_character_
    ),
    species   = factor(species),
    sal_adapt = factor(sal_adapt, levels = c(25,45), labels = c("A25","A45")),
    sal_exp   = factor(sal_exp,   levels = c(25,45), labels = c("E25","E45")),
    repl      = factor(repl)
  )

# Keep only DMS(P)t rows and t0/t24
df <- df0 %>% filter(id == "DMS(P)t", time_std %in% c("t0","t24"))

# --- Build paired dataset per bottle branch (control/shock at t24) ---
# Baseline at t0 (scaled to bottle start density)
base_t0 <- df %>%
  filter(time_std == "t0") %>%
  mutate(dens_t0_bottle = 0.6 * cell_dens_avr) %>%
  select(
    species, sal_adapt, repl,
    yield_t0 = yield,
    dens_t0  = dens_t0_bottle,
    size_t0  = cell_size_avr
  )

# t24 observed
t24 <- df %>%
  filter(time_std == "t24") %>%
  select(
    species, sal_adapt, sal_exp, repl,
    yield_t24 = yield,
    dens_t24  = cell_dens_avr,
    size_t24  = cell_size_avr
  )

# Join: each t24 branch gets its own t0 baseline
paired <- t24 %>%
  inner_join(base_t0, by = c("species","sal_adapt","repl")) %>%
  mutate(
    bottle    = interaction(species, sal_adapt, repl, drop = TRUE),
    d_yield   = yield_t24 - yield_t0,
    d_dens    = dens_t24  - dens_t0,
    d_size    = size_t24  - size_t0,
    pct_yield = ifelse(is.finite(yield_t0) & yield_t0 != 0,
                       (yield_t24 / yield_t0 - 1) * 100, NA_real_),
    pct_dens  = ifelse(is.finite(dens_t0)  & dens_t0  != 0,
                       (dens_t24  / dens_t0  - 1) * 100, NA_real_),
    pct_size  = ifelse(is.finite(size_t0)  & size_t0  != 0,
                       (size_t24  / size_t0  - 1) * 100, NA_real_)
  )

# --- Quick QC: how many paired branches per cell? (prints to console) ---
paired %>%
  count(species, sal_adapt, sal_exp, name = "n_branches") %>%
  arrange(species, sal_adapt, sal_exp) %>%
  print(n = 50)

# --- Analysis helper: Δ-LMM + emmeans (Holm) ---
analyze_delta <- function(dat, d_col, pct_col, nice_label, adjust = "holm") {
  d_use <- dat %>%
    select(species, sal_adapt, sal_exp, repl, bottle, !!d_col, !!pct_col) %>%
    rename(delta = !!d_col, pct = !!pct_col) %>%
    filter(!is.na(delta))
  if (nrow(d_use) < 8L) {
    message("[", nice_label, "] Too few paired rows; skipping."); return(NULL)
  }
  fit <- tryCatch(lmer(delta ~ species * sal_adapt * sal_exp + (1 | bottle), data = d_use),
                  error = function(e) NULL)
  if (is.null(fit)) {
    message("[", nice_label, "] lmer failed; trying lm without random effect.")
    fit <- tryCatch(lm(delta ~ species * sal_adapt * sal_exp, data = d_use),
                    error = function(e) NULL)
  }
  if (is.null(fit)) { message("[", nice_label, "] Modeling failed."); return(NULL) }
  
  aov_tbl <- car::Anova(fit, type = 3) %>% broom::tidy()
  
  out_sets <- list(
    "Species within sal_adapt × sal_exp" =
      pairs(emmeans(fit, ~ species   | sal_adapt * sal_exp), adjust = adjust),
    "E25 vs E45 within species × sal_adapt" =
      pairs(emmeans(fit, ~ sal_exp   | species * sal_adapt), adjust = adjust),
    "Acclimation (A25 vs A45) within species × sal_exp" =
      pairs(emmeans(fit, ~ sal_adapt | species * sal_exp),   adjust = adjust),
    "Species (marginal)"   = pairs(emmeans(fit, ~ species),   adjust = adjust),
    "Acclimation (marginal)" = pairs(emmeans(fit, ~ sal_adapt), adjust = adjust),
    "Exposure (marginal)"     = pairs(emmeans(fit, ~ sal_exp),   adjust = adjust),
    "All 3-way cells (Tukey)" = pairs(emmeans(fit, ~ species * sal_adapt * sal_exp),
                                      adjust = "tukey")
  )
  
  contrasts_df <- lapply(names(out_sets), function(ttl) {
    as.data.frame(summary(out_sets[[ttl]], infer = c(TRUE, TRUE))) %>%
      mutate(contrast_set = ttl, .before = 1)
  }) %>%
    bind_rows() %>%
    mutate(label = nice_label,
           id = "Physiology Δ", response = as.character(d_col))
  
  cell_means <- as.data.frame(summary(emmeans(fit, ~ species * sal_adapt * sal_exp), infer = TRUE)) %>%
    mutate(label = nice_label, id = "Physiology Δ", response = as.character(d_col))
  
  pct_tbl <- d_use %>%
    filter(!is.na(pct)) %>%
    group_by(species, sal_adapt, sal_exp) %>%
    summarise(
      mean_pct = mean(pct, na.rm = TRUE),
      se_pct   = sd(pct, na.rm = TRUE) / sqrt(sum(!is.na(pct))),
      n        = n(), .groups = "drop"
    ) %>%
    mutate(label = nice_label, id = "Physiology Δ", response = paste0(as.character(d_col), "_pct"))
  
  list(fit = fit, anova = aov_tbl, contrasts = contrasts_df,
       cell_means = cell_means, pct = pct_tbl, label = nice_label)
}

# --- Run analyses for 3 metrics ---
res_size  <- analyze_delta(paired, sym("d_size"),  sym("pct_size"),  "Physiology Δ: Mean cell size")
res_dens  <- analyze_delta(paired, sym("d_dens"),  sym("pct_dens"),  "Physiology Δ: Cell density")
res_yield <- analyze_delta(paired, sym("d_yield"), sym("pct_yield"), "Physiology Δ: Yield (Fv/Fm)")

# Drop NULLs and name them
physio_results <- Filter(Negate(is.null), list(
  physio_delta_size    = res_size,
  physio_delta_density = res_dens,
  physio_delta_yield   = res_yield
))

# --- Merge into your main 'results' list (create if it doesn't exist) ---
if (!exists("results")) results <- list()
results <- c(results, physio_results)

# Keep only successful fits with a non-empty contrasts table
ok_names <- names(results)[vapply(
  results,
  function(x) !is.null(x) && !is.null(x$contrasts) && nrow(x$contrasts) > 0,
  logical(1)
)]
stopifnot(length(ok_names) > 0)

# --- Export to the SAME combined files as the rest ---
out_dir <- get0("out_dir", ifnotfound = "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
combined_base <- file.path(out_dir, "ALL_contrasts_combined")

make_family_table_tex_all(ok_names, results, combined_base)                 # LaTeX (all contrasts)
make_family_txt_all(ok_names, results, paste0(combined_base, ".txt"))       # TXT   (all contrasts)

# --- Optional quick screen plot (Δ by cell) ---
if (nrow(paired)) {
  long_d <- paired %>%
    select(species, sal_adapt, sal_exp, repl, d_yield, d_dens, d_size) %>%
    pivot_longer(starts_with("d_"), names_to = "metric", values_to = "delta") %>%
    mutate(metric = dplyr::recode(
      metric,
      d_yield = "Yield (Δ)", d_dens = "Cell density (Δ)", d_size = "Mean cell size (Δ)"
    ))
  p <- ggplot(long_d, aes(sal_exp, delta)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.6, size = 1.7) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    facet_grid(metric ~ species + sal_adapt, scales = "free_y") +
    labs(x = "Exposure (t24)", y = "Δ(t24 − t0)", title = "Physiology changes by cell") +
    theme_minimal(base_size = 12)
  print(p)
  # ggsave(file.path(out_dir, "physiology_t0t24_DMSPt_QC.png"), p, width = 8, height = 6, dpi = 300)
}

# =====================================================
# BOX PLOTS — t24 (standard style) for all metrics
#   - Uses your ALT style for uptake_alt & loss_alt
#   - One PNG per metric + one multipage PDF
#   - Optional Holm CLD letters (off by default)
# =====================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(purrr); library(rlang)
  library(emmeans)   # only used if add_letters = TRUE
})

# ---- Output dir for figures ----
fig_dir <- file.path(getwd(), "results", "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+$", "", x); tolower(x)
}
setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")

# ---- Load & harmonise once (same as your main script) ----
chr <- read_csv("clean_chr.csv",           show_col_types = FALSE)
pol <- read_csv("clean_pol_modified.csv",  show_col_types = FALSE)

dat_all <- bind_rows(chr, pol) %>%
  mutate(
    time_std = case_when(
      as.character(time) %in% c("0","t0")   ~ "t0",
      as.character(time) %in% c("24","t24") ~ "t24",
      TRUE ~ NA_character_
    ),
    species   = factor(species),
    sal_adapt = factor(sal_adapt, levels = c(25,45), labels = c("A25","A45")),
    sal_exp   = factor(sal_exp,   levels = c(25,45), labels = c("E25","E45")),
    repl      = factor(repl)
  )

# =====================================================
# Build ALT metrics (as in your example)
#   uptake_alt = (d3_p_taken_up × μ_POC) / POC_t24
#   loss_alt   = (d3_p_lost_demethylated × μ_POC) / POC_t24
# =====================================================
poc_mu_t24 <- dat_all %>%
  filter(id == "DMSPp", time_std == "t24") %>%
  select(species, sal_adapt, sal_exp, repl,
         POC_t24 = total_poc_mg_l,
         mu_POC  = u_poc_h_1)

upt_loss_t24 <- dat_all %>%
  filter(id == "DMS(P)t", time_std == "t24") %>%
  select(species, sal_adapt, sal_exp, repl,
         uptake = d3_p_taken_up,
         loss   = d3_p_lost_demethylated)

gc_alt <- upt_loss_t24 %>%
  left_join(poc_mu_t24, by = c("species","sal_adapt","sal_exp","repl")) %>%
  mutate(
    uptake_alt = if_else(is.finite(mu_POC) & is.finite(POC_t24) & POC_t24 > 0,
                         (uptake * mu_POC) / POC_t24, NA_real_),
    loss_alt   = if_else(is.finite(mu_POC) & is.finite(POC_t24) & POC_t24 > 0,
                         (loss   * mu_POC) / POC_t24, NA_real_)
  )

# =====================================================
# Plot helper (your standard style) + optional CLD letters
# =====================================================
make_box <- function(dat, value_col, y_label, add_letters = FALSE) {
  dat_p <- dat %>%
    mutate(
      sal_adapt = factor(sal_adapt),
      sal_exp   = factor(sal_exp),
      group_x   = interaction(species, sal_adapt, sep = " × ")
    )
  p <- ggplot(dat_p, aes(x = group_x, y = .data[[value_col]], fill = sal_exp)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(aes(color = sal_exp), width = 0.15, alpha = 0.6,
                size = 1.6, show.legend = FALSE) +
    labs(x = "Species × acclimation salinity", y = y_label, fill = "Exposure") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),
          legend.position = "top")
  
  if (isTRUE(add_letters)) {
    try({
      mod <- lm(as.formula(paste(value_col, "~ species * sal_adapt * sal_exp")),
                data = dat_p)
      emm <- emmeans(mod, ~ species * sal_adapt * sal_exp)
      cld_df <- emmeans::cld(emm, adjust = "holm", Letters = letters) %>%
        as.data.frame() %>%
        mutate(group_x = interaction(species, sal_adapt, sep = " × "))
      y_pos <- dat_p %>%
        group_by(group_x) %>%
        summarise(y_top = max(.data[[value_col]], na.rm = TRUE) * 1.05,
                  .groups = "drop")
      cld_df <- left_join(cld_df, y_pos, by = "group_x")
      p <- p + geom_text(data = cld_df,
                         aes(x = group_x, y = y_top, label = .group),
                         inherit.aes = FALSE)
    }, silent = TRUE)
  }
  p
}


# =====================================================
# SPEC: what to plot (direct-from-column metrics @ t24)
# id = source pool; col = column name; label = axis label
# Any missing columns are skipped with a message.
# =====================================================
direct_specs <- tribble(
  ~id,       ~col,                        ~label,
  "DMSPp",   "u_poc_h_1",                 "μ-POC (h⁻¹)",
  "DMSPp",   "u_dmsp_h_1",                "μ-DMSP (h⁻¹) from DMSPp",
  "DMS",     "u_dmsp_h_1",                "μ-DMSP (h⁻¹) from DMS",
  "DMSPp",   "u_dmsp_u_poc",              "μ-DMSP / μ-POC",
  "DMSPp",   "dmsp_c_poc_mol_mol",        "DMSP-C : POC (mol:mol)",
  "DMS(P)t", "d3_p_taken_up",             "D3-DMSP uptake (raw)",
  "DMS(P)t", "d3_p_lost_demethylated",    "D3-DMSP loss (demethyl.)",
  "DMS(P)t", "dms_puptake_c_poc_mol_mol", "DMSP uptake-C : POC (mol:mol)",
  "DMS(P)t", "dmsp_uptake_of_total_dmsp", "DMSP uptake (% of total DMSP)"  # will be skipped if absent
  # If you want: "DMSPp","up_dmspp","up_dmspp (units?)"
)

# ---- Generate & save: direct metrics ----
plots <- list()
for (i in seq_len(nrow(direct_specs))) {
  spec <- direct_specs[i,]
  if (!spec$col %in% names(dat_all)) {
    message("Skipping '", spec$col, "' (column not found)."); next
  }
  df_plot <- dat_all %>%
    filter(id == spec$id, time_std == "t24") %>%
    select(species, sal_adapt, sal_exp, repl, value = !!sym(spec$col))
  if (!nrow(df_plot)) { message("No t24 rows for ", spec$col, " in ", spec$id); next }
  p <- make_box(df_plot, "value", spec$label, add_letters = FALSE)
  fname <- file.path(fig_dir, paste0("box_", safe_filename(paste(spec$id, spec$col, sep = "_")), ".png"))
  ggsave(fname, p, width = 7, height = 4.5, dpi = 300, bg = "white")
  plots[[paste(spec$id, spec$col, sep = "_")]] <- p
}

# ---- ALT metrics (from gc_alt built above) ----
alt_specs <- tribble(
  ~col,         ~label,                                      ~stub,
  "uptake_alt", "(D3P-Uptake × μ-POC) / POC(t24)",          "alt_uptake",
  "loss_alt",   "(D3P-Loss × μ-POC) / POC(t24)",            "alt_loss"
)

for (i in seq_len(nrow(alt_specs))) {
  spec <- alt_specs[i,]
  if (!spec$col %in% names(gc_alt)) next
  df_plot <- gc_alt %>% select(species, sal_adapt, sal_exp, repl, value = !!sym(spec$col))
  p <- make_box(df_plot, "value", spec$label, add_letters = FALSE)
  fname <- file.path(fig_dir, paste0("box_", spec$stub, ".png"))
  ggsave(fname, p, width = 7, height = 4.5, dpi = 300, bg = "white")
  plots[[spec$stub]] <- p
}

# ---- Also write a single multipage PDF with all plots ----
if (length(plots)) {
  pdf(file.path(fig_dir, "all_boxplots_t24.pdf"), width = 7, height = 4.5, onefile = TRUE)
  for (nm in names(plots)) print(plots[[nm]])
  dev.off()
  message("Saved: ", normalizePath(file.path(fig_dir, "all_boxplots_t24.pdf"), winslash = "/"))
} else {
  message("No plots produced.")
}
