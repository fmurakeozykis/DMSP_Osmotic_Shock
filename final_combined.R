# ============================================================
# ONE-BLOCK PIPELINE (prose TXT only; no LaTeX)
# Frame = first script (skew/log1p, AIC selection incl. lme+varIdent)
# Metrics = second script (incl. normalized_df)
# Inference = Type-III on clean lmer refit; if lmer is SINGULAR → drop RE and use lm
# Post-hoc = GATED (only if relevant omnibus terms are significant)
# Output = ONE TXT with transparent prose (assumptions + changes)# ============================================================
# ONE-BLOCK PIPELINE (prose TXT only; no LaTeX)
# t24-only changes: use raw columns directly by ID (no re-compute of responses)
# ============================================================
set.seed(42)

quick_diag <- function(fit, data = NULL, grp = NULL, label = NULL) {
  is_lmer <- inherits(fit, "merMod")
  res <- tryCatch(residuals(fit, type = if (is_lmer) "pearson" else "response"),
                  error = function(e) residuals(fit))
  fv  <- tryCatch(fitted(fit), error = function(e) rep(NA_real_, length(res)))
  n   <- sum(is.finite(res))
  idx <- which(is.finite(res)); idx <- if (length(idx) > 5000) sample(idx, 5000) else idx
  shap_p <- if (length(idx) >= 3) tryCatch(shapiro.test(res[idx])$p.value, error=function(e) NA_real_) else NA_real_
  sk   <- tryCatch(e1071::skewness(res, na.rm=TRUE), error=function(e) NA_real_)
  kur  <- tryCatch(e1071::kurtosis(res, na.rm=TRUE), error=function(e) NA_real_)
  lev_p <- NA_real_
  if (!is.null(data) && all(c("species","sal_adapt","sal_exp") %in% names(data))) {
    grp_fac <- with(data, interaction(species, sal_adapt, sal_exp, drop=TRUE))
    lev_p <- tryCatch(car::leveneTest(res ~ grp_fac, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
  }
  ncv_p <- if (!is_lmer) tryCatch(car::ncvTest(fit)$`p`, error=function(e) NA_real_) else NA_real_
  sp_p  <- tryCatch(cor.test(abs(res), fv, method="spearman")$p.value, error=function(e) NA_real_)
  n_std3 <- tryCatch(sum(abs(scale(res)) > 3, na.rm=TRUE), error=function(e) NA_integer_)
  max_cook <- NA_real_; n_cook <- NA_integer_
  if (!is_lmer) {
    cd <- tryCatch(cooks.distance(fit), error=function(e) rep(NA_real_, n))
    thr <- 4 / n
    max_cook <- suppressWarnings(max(cd, na.rm=TRUE))
    n_cook   <- sum(cd > thr, na.rm=TRUE)
  }
  re_var <- NA_real_; resid_var <- NA_real_; icc <- NA_real_; singular <- NA
  if (is_lmer) {
    vc <- tryCatch(lme4::VarCorr(fit), error=function(e) NULL)
    if (!is.null(vc)) {
      resid_var <- tryCatch(attr(vc, "sc")^2, error=function(e) NA_real_)
    }
    singular <- tryCatch(lme4::isSingular(fit, tol=1e-4), error=function(e) NA)
  }
  lines <- c(
    if (!is.null(label)) paste0("Diagnostics for: ", label),
    sprintf("  n(resid)=%s; Shapiro p=%s; skew=%.2f; kurt=%.2f",
            n, ifelse(is.na(shap_p), "NA", formatC(shap_p, digits=3)),
            ifelse(is.na(sk), NaN, sk), ifelse(is.na(kur), NaN, kur)),
    sprintf("  Heteroskedasticity — Levene p=%s; %s; Spearman(|res|,fitted) p=%s",
            ifelse(is.na(lev_p), "NA", formatC(lev_p, digits=3)),
            if (!is_lmer) paste0("ncvTest(BP) p=", ifelse(is.na(ncv_p), "NA", formatC(ncv_p, digits=3))) else "ncvTest(BP) (lm only)",
            ifelse(is.na(sp_p), "NA", formatC(sp_p, digits=3))),
    sprintf("  Outliers — |std resid|>3: %s; %s",
            ifelse(is.na(n_std3), "NA", n_std3),
            if (!is_lmer) paste0("Cook's: max=", formatC(max_cook, digits=3),
                                 ", #>4/n=", ifelse(is.na(n_cook), "NA", n_cook))
            else "Cook's (lm only)")
  )
  list(lines = lines)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
  library(lme4); library(nlme); library(emmeans); library(car); library(MuMIn)
  library(broom); library(e1071)
})

# ---------- options ----------
set.seed(1)
.old_contr <- options("contrasts")
options(contrasts = c("contr.sum","contr.poly"))
on.exit(do.call(options, .old_contr), add = TRUE)
ALLOW_SINGULAR <- TRUE
alpha <- 0.05
`%||%` <- function(a,b) if(!is.null(a)) a else b
fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 1e-3, formatC(p, format="e", digits=2), sprintf("%.3f", p)))
fmt_n <- function(x, d=3) ifelse(is.na(x), "NA", formatC(x, format="f", digits=d))
hdr <- function(txt){ c("", paste0("## ", txt), strrep("-", nchar(txt)+3)) }

# ---------- load & harmonise (adjust paths if needed) ----------
chr <- read_csv("clean_chr.csv", show_col_types = FALSE)
pol <- tryCatch(read_csv("clean_pol_modified.csv", show_col_types = FALSE),
                error = function(e) read_csv("clean_pol.csv", show_col_types = FALSE))

dat_all <- bind_rows(chr, pol) %>%
  mutate(
    time_num = suppressWarnings(as.integer(as.character(time))),
    time_std = if_else(time_num == 0L, "t0", "t24"),
    time     = factor(time_std, levels = c("t0","t24")),
    species   = factor(species),
    sal_adapt = factor(sal_adapt, levels = c(25,45), labels = c("A25","A45")),
    sal_exp   = factor(sal_exp,   levels = c(25,45), labels = c("E25","E45")),
    repl      = factor(repl),
    bottle    = interaction(species, sal_adapt, repl, drop = TRUE),
    id        = factor(id)
  ) %>%
  select(-time_num)

# ---------- t24-only normalized bundle (direct pulls; no re-compute of core responses) ----------
keys <- c("species","sal_adapt","sal_exp","repl")

# DMSPp rows (t24)
p24 <- dat_all %>%
  filter(time_std == "t24", id == "DMSPp") %>%
  select(
    all_of(keys),
    mu_POC = u_poc_h_1,               # μ-POC (direct)
    u_dmsp_u_poc,                     # μ-DMSP/μ-POC (direct)
    dmsp_c_poc_mol_mol,               # DMSP-C:POC (direct)
    any_of(c(
      "up_dmspp",                     # Uptake % of DMSPp (direct)
      "dmsp_uptake_of_total_dmsp",    # Fraction of total DMSP (direct)
      "dms_puptake_c_poc_mol_mol"     # Uptake-C:POC (direct)
    ))
  )

# DMS(P)t rows (t24): direct tracer uptake/loss
t24 <- dat_all %>%
  filter(time_std == "t24", id == "DMS(P)t") %>%
  select(
    all_of(keys),
    uptake_alt = d3_p_taken_up,           # direct
    loss_alt   = d3_p_lost_demethylated   # direct
  )

# DMSPd rows (t24): direct μ-DMSP-from-DMSPd (used as-is and/or for ratios if you want them elsewhere)
d24 <- dat_all %>%
  filter(time_std == "t24", id == "DMSPd") %>%
  select(all_of(keys), mu_DMSPd = u_dmsp_h_1)   # direct

# Assemble a t24 "normalized" table that only *joins* direct columns
normalized_df <- p24 %>%
  left_join(t24, by = keys) %>%
  left_join(d24, by = keys) %>%
  mutate(
    time_std = "t24",
    time     = factor("t24", levels = c("t0","t24")),
    id       = "Normalized"
  )

# ---------- analysis list (each pulls t24 + correct ID; responses taken *directly*) ----------
analyses <- list(
  # Production metrics (by source id)
  list(name="mu_POC_DMSPp",           dat=dat_all, y="u_poc_h_1",              id="DMSPp",   label="μ-POC (from DMSPp)"),
  list(name="mu_DMSP_from_DMSPp",     dat=dat_all, y="u_dmsp_h_1",             id="DMSPp",   label="μ-DMSP (from DMSPp)"),
  list(name="mu_DMSP_from_DMS",       dat=dat_all, y="u_dmsp_h_1",             id="DMS",     label="μ-DMSP (from DMS)"),
  list(name="mu_DMSP_from_DMSPd",     dat=dat_all, y="u_dmsp_h_1",             id="DMSPd",   label="μ-DMSP (from DMSPd)"),
  
  # Tracer uptake/loss (direct)
  list(name="uptake_d3",              dat=dat_all, y="d3_p_taken_up",          id="DMS(P)t", label="D3-DMSP taken up"),
  list(name="loss_d3",                dat=dat_all, y="d3_p_lost_demethylated", id="DMS(P)t", label="D3-DMSP lost (demethylated)"),
  
  # Normalized table (still direct pulls, just joined)
  list(name="norm_u_dmsp_u_poc",      dat=normalized_df, y="u_dmsp_u_poc",           id="Normalized", label="μ-DMSP / μ-POC (direct)"),
  list(name="norm_dmsp_c_poc",        dat=normalized_df, y="dmsp_c_poc_mol_mol",     id="Normalized", label="DMSP-C : POC (mol:mol)"),
  list(name="norm_uptake_alt_direct", dat=normalized_df, y="uptake_alt",             id="Normalized", label="Tracer uptake (d3_p_taken_up)"),
  list(name="norm_loss_alt_direct",   dat=normalized_df, y="loss_alt",               id="Normalized", label="Tracer loss (d3_p_lost_demethylated)"),
  list(name="norm_mu_DMSPd_direct",   dat=normalized_df, y="mu_DMSPd",               id="Normalized", label="μ-DMSP (from DMSPd; direct)"),
  
  # Exact extras on DMSPp ID (t24)
  list(name="up_dmspp_t24",           dat=dat_all, y="up_dmspp",                  id="DMSPp", label="Uptake (% of DMSPp)"),
  list(name="dmsp_upt_total_t24",     dat=dat_all, y="dmsp_uptake_of_total_dmsp", id="DMSPp", label="DMSP uptake of total DMSP"),
  list(name="dms_puptake_c_poc_t24",  dat=dat_all, y="dms_puptake_c_poc_mol_mol", id="DMSPp", label="DMSP uptake-C : POC (mol:mol)")
)

# ---------- core engines ----------
refit_lmer <- function(data, resp){
  form <- as.formula(paste(resp, "~ species * sal_adapt * sal_exp + (1|bottle)"))
  lme4::lmer(form, data = data, REML = TRUE,
             control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

analyze_one <- function(dat, y, label, id_filter, time_filter = "t24"){
  lines <- c(hdr(label))
  
  d0 <- dat %>%
    filter(time_std == time_filter, id == id_filter) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp, repl) %>%
    mutate(bottle = interaction(species, sal_adapt, repl, drop = TRUE)) %>%
    mutate(.y = suppressWarnings(as.numeric(.data[[y]]))) %>%
    filter(is.finite(.y))
  
  if (nrow(d0) < 8L ||
      n_distinct(d0$species)   < 2 ||
      n_distinct(d0$sal_adapt) < 2 ||
      n_distinct(d0$sal_exp)   < 2) {
    return(c(lines, "Skipped: insufficient data or factors lack ≥2 levels.\n"))
  }
  if (length(unique(d0$.y)) < 2) {
    return(c(lines, "Skipped: response has < 2 unique finite values.\n"))
  }
  
  base_fit <- tryCatch(refit_lmer(d0, ".y"), error=function(e) NULL)
  if (is.null(base_fit)) return(c(lines, "Skipped: baseline lmer failed.\n"))
  
  r <- residuals(base_fit, type="pearson")
  shapiro_p <- tryCatch(shapiro.test(as.numeric(r))$p.value, error=function(e) NA_real_)
  grp <- with(d0, interaction(species, sal_adapt, sal_exp, drop = TRUE))
  levene_p <- tryCatch(car::leveneTest(r ~ grp, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
  
  assumptions <- c(
    paste0("Normality (Shapiro–Wilk on residuals): p = ", fmt_p(shapiro_p)),
    paste0("Homoskedasticity (Levene across S×A×E): p = ", fmt_p(levene_p))
  )
  flagged <- c()
  if (!is.na(shapiro_p) && shapiro_p < alpha) flagged <- c(flagged, paste0("Normality (p=", fmt_p(shapiro_p), ")"))
  if (!is.na(levene_p)  && levene_p  < alpha) flagged <- c(flagged,  paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
  
  changes <- character()
  transform <- NULL
  y_use <- ".y"
  
  # Optional transform only if strongly skewed (no re-compute of response itself)
  sk <- suppressWarnings(e1071::skewness(d0$.y, na.rm = TRUE))
  if (is.finite(sk) && sk > 1 && all(d0$.y >= 0, na.rm = TRUE)) {
    d0$.y_log1p <- log1p(d0$.y)
    base_fit2 <- tryCatch(refit_lmer(d0, ".y_log1p"), error=function(e) NULL)
    if (!is.null(base_fit2)) {
      base_fit <- base_fit2; transform <- "log1p"; y_use <- ".y_log1p"
      changes <- c(changes, "Applied log1p due to skew > 1; post-hoc back-transformed.")
      r2 <- residuals(base_fit, type="pearson")
      levene_p <- tryCatch(car::leveneTest(r2 ~ grp, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
      flagged <- flagged[!grepl("^Homoskedasticity", flagged)]
      if (!is.na(levene_p) && levene_p < alpha)
        flagged <- c(flagged, paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
    } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
      changes <- c(changes, "Normality flagged but log1p not applied (failed refit or negative values).")
    }
  } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
    changes <- c(changes, "Normality flagged; retained scale (no positive-only transform justified).")
  }
  
  # Selection block for heteroskedasticity if needed
  selected <- base_fit
  if (!is.na(levene_p) && levene_p < alpha) {
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * sal_exp"))
    fit_try <- function(w) tryCatch(nlme::lme(fixed=form_fixed, random=~1|bottle, data=d0,
                                              method="REML", weights=w, control=nlme::lmeControl(msMaxIter=200)),
                                    error=function(e) NULL)
    cands <- list(
      none = fit_try(NULL),
      byE  = fit_try(nlme::varIdent(~1|sal_exp)),
      byS  = fit_try(nlme::varIdent(~1|species)),
      byA  = fit_try(nlme::varIdent(~1|sal_adapt)),
      byAE = fit_try(nlme::varIdent(~1|sal_adapt:sal_exp))
    )
    ok <- cands[!vapply(cands, is.null, logical(1))]
    if (length(ok)) {
      aics <- sapply(ok, AIC); pick <- names(which.min(aics))
      selected <- ok[[pick]]
      changes <- c(changes, paste0("Modeled heteroskedasticity via nlme::lme + varIdent (", pick, ") by AIC."))
    }
  }
  
  # Inference engine: lmer unless singular → lm
  fit_inf <- tryCatch(refit_lmer(d0, y_use), error=function(e) NULL)
  lines <- c(lines, "Diagnostics (auto, no plots):",
             quick_diag(fit_inf, data = d0, label = paste0(label, " (", id_filter, " @ ", time_filter, ")"))$lines)
  
  if (is.null(fit_inf)) return(c(lines, "Skipped: inference refit failed.\n"))
  used_engine <- "lmer"
  if (inherits(fit_inf, "merMod") && lme4::isSingular(fit_inf, tol=1e-4)) {
    fit_inf <- stats::lm(as.formula(paste(y_use, "~ species * sal_adapt * sal_exp")), data = d0)
    used_engine <- "lm"
    changes <- c(changes, "Random intercept singular; dropped RE and used OLS (lm) for inference.")
  }
  
  aov_tab <- suppressWarnings(car::Anova(fit_inf, type="III"))
  aov_df  <- as.data.frame(aov_tab)
  pcol <- intersect(colnames(aov_df), c("Pr(>F)","Pr(>Chisq)","Pr(>Chi)","Pr..F."))
  if (!length(pcol)) return(c(lines, "Type-III ANOVA failed (no p-value column).", ""))
  
  getp <- function(term){
    row <- aov_df[rownames(aov_df)==term, , drop=FALSE]
    if (!nrow(row)) return(NA_real_)
    as.numeric(row[[pcol[1]]][1])
  }
  p3  <- getp("species:sal_adapt:sal_exp")
  pSA <- getp("species:sal_adapt")
  pSE <- getp("species:sal_exp")
  pAE <- getp("sal_adapt:sal_exp")
  pS  <- getp("species")
  pA  <- getp("sal_adapt")
  pE  <- getp("sal_exp")
  
  R2 <- suppressWarnings(tryCatch(MuMIn::r.squaredGLMM(selected), error=function(e) c(NA_real_, NA_real_)))
  
  lines <- c(lines,
             "Assumption checks:",
             paste0("  - ", assumptions),
             paste0("Assumptions flagged (α=", alpha, "): ",
                    if (length(flagged)) paste(flagged, collapse = "; ") else "none."),
             if (length(changes)) c("Changes made:", paste0("  - ", changes)) else "Changes made: none.",
             paste0("Inference engine: ", if (used_engine=="lm") "OLS (lm) — random intercept dropped" else "lmer (random-intercept)"),
             "Type-III ANOVA (on inference engine above):",
             paste0("  - species: p=", fmt_p(pS),
                    "; sal_adapt: p=", fmt_p(pA),
                    "; sal_exp: p=", fmt_p(pE)),
             paste0("  - species:sal_adapt: p=", fmt_p(pSA),
                    "; species:sal_exp: p=", fmt_p(pSE),
                    "; sal_adapt:sal_exp: p=", fmt_p(pAE)),
             paste0("  - species:sal_adapt:sal_exp: p=", fmt_p(p3)),
             paste0("Model R² (selected fit) — marginal=", fmt_n(as.numeric(R2[1])),
                    ", conditional=", fmt_n(as.numeric(R2[2]))))
  
  proceed <- any(c(p3,pSA,pSE,pAE,pS,pA,pE) < alpha, na.rm = TRUE)
  if (!proceed) return(c(lines, "Post-hoc: skipped (no omnibus terms significant at α=0.05).", ""))
  
  trn <- NULL
  regrid_if <- function(em) if (!is.null(trn) && trn=="log1p") emmeans::regrid(em, transform="response") else em
  posthoc_lines <- character()
  
  if (any(c(p3,pSA,pSE,pAE) < alpha, na.rm = TRUE)) {
    em1 <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_exp, by = c("species","sal_adapt")))
    cmp1 <- summary(emmeans::contrast(em1, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    as_lines <- function(df) if (!nrow(df)) character() else
      sprintf("  - %s: %s [%s, %s]; p=%s",
              df$contrast, fmt_n(df$estimate), fmt_n(df$lower.CL), fmt_n(df$upper.CL), fmt_p(df$p.value))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (E within species × acclimation):", as_lines(as.data.frame(cmp1)))
    
    emS <- regrid_if(emmeans::emmeans(fit_inf, ~ species, by = c("sal_adapt","sal_exp")))
    cmpS <- summary(emmeans::contrast(emS, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Species within acclimation × exposure):", as_lines(as.data.frame(cmpS)))
    
    emA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt, by = c("species","sal_exp")))
    cmpA <- summary(emmeans::contrast(emA, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Acclimation within species × exposure):", as_lines(as.data.frame(cmpA)))
    
    em_all <- regrid_if(emmeans::emmeans(fit_inf, ~ species * sal_adapt * sal_exp))
    cmp_all <- summary(pairs(em_all, adjust="tukey"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (All 3-way cell means; Tukey):", as_lines(as.data.frame(cmp_all)))
  } else {
    if (!is.na(pS) && pS < alpha) {
      em_mS <- regrid_if(emmeans::emmeans(fit_inf, ~ species))
      cmp_mS <- summary(pairs(em_mS, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Species; marginal):",
                         sprintf("  - %s: %s [%s, %s]; p=%s",
                                 cmp_mS$contrast, fmt_n(cmp_mS$estimate), fmt_n(cmp_mS$lower.CL),
                                 fmt_n(cmp_mS$upper.CL), fmt_p(cmp_mS$p.value)))
    }
    if (!is.na(pA) && pA < alpha) {
      em_mA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt))
      cmp_mA <- summary(pairs(em_mA, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Acclimation; marginal):",
                         sprintf("  - %s: %s [%s, %s]; p=%s",
                                 cmp_mA$contrast, fmt_n(cmp_mA$estimate), fmt_n(cmp_mA$lower.CL),
                                 fmt_n(cmp_mA$upper.CL), fmt_p(cmp_mA$p.value)))
    }
    if (!is.na(pE) && pE < alpha) {
      em_mE <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_exp))
      cmp_mE <- summary(pairs(em_mE, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Exposure; marginal):",
                         sprintf("  - %s: %s [%s, %s]; p=%s",
                                 cmp_mE$contrast, fmt_n(cmp_mE$estimate), fmt_n(cmp_mE$lower.CL),
                                 fmt_n(cmp_mE$upper.CL), fmt_p(cmp_mE$p.value)))
    }
  }
  c(lines, if (length(posthoc_lines)) posthoc_lines else "Post-hoc: none produced.", "")
}

# ---------- Run all analyses & write TXT ----------
prose <- c("==== MIXED-MODEL RESULTS (prose) ====",
           paste0("Alpha = ", alpha),
           "Type-III ANOVA computed on a clean lmer refit; if the random intercept is singular,",
           "we drop it and use OLS (lm) for inference on fixed effects. Post-hoc uses the same engine.",
           "t24-only: responses pulled directly from requested columns by ID; no re-computation of μ metrics.",
           "")

for (a in analyses) {
  if (!a$y %in% names(a$dat)) {
    prose <- c(prose, hdr(a$label), "Skipped: response column not present.\n")
    next
  }
  blk <- tryCatch(analyze_one(a$dat, a$y, a$label, a$id, time_filter = "t24"),
                  error=function(e) c(hdr(a$label), paste("Analysis failed:", conditionMessage(e)), ""))
  prose <- c(prose, blk)
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path("results", "ANOVA_prose_results.txt")
writeLines(prose, out_txt)
cat("Wrote TXT to: ", normalizePath(out_txt, winslash = "/"), "\n", sep="")

# (Physiology Δ-pipeline remains as in your original script; not modified here.)




# ============================================================
# PHYSIOLOGY (t0 vs t24; exposure only at t24, NoExp is reference)
# Key change: E_phase is treatment-coded to avoid aliasing with time.
# ============================================================

phys_metrics <- c("cell_size_avr","yield","cell_dens_avr")
phys_labels  <- c(cell_size_avr="Cell size (avg)",
                  yield        ="Yield",
                  cell_dens_avr="Cell density (avg)")

# Build physiology dataframe (no exposure at t0; do NOT match acclimation)
d_phys <- dat_all %>%
  dplyr::filter(id == "DMS(P)t") %>%
  mutate(
    time   = factor(time, levels = c("t0","t24")),
    E_phase = factor(ifelse(time == "t24", as.character(sal_exp), "NoExp"),
                     levels = c("NoExp","E25","E45")),
    bottle = interaction(species, sal_adapt, repl, drop = TRUE)
  )

# IMPORTANT: treatment contrasts ONLY for E_phase (NoExp = reference)
d_phys$E_phase_trt <- d_phys$E_phase
contrasts(d_phys$E_phase_trt) <- stats::contr.treatment(nlevels(d_phys$E_phase_trt))

refit_lmer_phys <- function(data, resp){
  form <- as.formula(paste(resp, "~ species * sal_adapt * time + time:E_phase_trt + (1|bottle)"))
  lme4::lmer(form, data = data, REML = TRUE,
             control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

analyze_phys_one <- function(dat, y, label){
  lines <- c(hdr(label),
             "Physiology metrics are analyzed on raw columns with a time factor (t0 & t24).",
             "Exposure is coded only at t24 (E_phase = NoExp/E25/E45; NoExp is reference).")
  
  d0 <- dat %>%
    tidyr::drop_na(species, sal_adapt, time) %>%
    mutate(.y = suppressWarnings(as.numeric(.data[[y]]))) %>%
    filter(is.finite(.y))
  
  if (nrow(d0) < 8L ||
      dplyr::n_distinct(d0$species)   < 2 ||
      dplyr::n_distinct(d0$sal_adapt) < 2 ||
      dplyr::n_distinct(d0$time)      < 2) {
    return(c(lines, "Skipped: insufficient data or a factor lacks ≥2 levels.\n"))
  }
  if (length(unique(d0$.y)) < 2) {
    return(c(lines, "Skipped: response has < 2 unique finite values.\n"))
  }
  
  base_fit <- tryCatch(refit_lmer_phys(d0, ".y"), error=function(e) NULL)
  if (is.null(base_fit)) return(c(lines, "Skipped: baseline lmer failed.\n"))
  
  # Diagnostics — Levene across S×A×E_phase (original coding for grouping)
  r <- residuals(base_fit, type="pearson")
  shapiro_p <- tryCatch(shapiro.test(as.numeric(r))$p.value, error=function(e) NA_real_)
  grp3 <- with(d0, interaction(species, sal_adapt, E_phase, drop = TRUE))
  levene_p <- tryCatch(car::leveneTest(r ~ grp3, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
  
  assumptions <- c(
    paste0("  - Normality (Shapiro–Wilk on residuals): p = ", fmt_p(shapiro_p)),
    paste0("  - Homoskedasticity (Levene across S×A×E_phase): p = ", fmt_p(levene_p))
  )
  flagged <- c()
  if (!is.na(shapiro_p) && shapiro_p < alpha) flagged <- c(flagged, paste0("Normality (p=", fmt_p(shapiro_p), ")"))
  if (!is.na(levene_p)  && levene_p  < alpha) flagged <- c(flagged,  paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
  
  changes <- character()
  transform <- NULL
  y_use <- ".y"
  
  sk <- suppressWarnings(e1071::skewness(d0$.y, na.rm = TRUE))
  if (is.finite(sk) && sk > 1 && all(d0$.y >= 0, na.rm = TRUE)) {
    d0$.y_log1p <- log1p(d0$.y)
    base_fit2 <- tryCatch(refit_lmer_phys(d0, ".y_log1p"), error=function(e) NULL)
    if (!is.null(base_fit2)) { base_fit <- base_fit2; transform <- "log1p"; y_use <- ".y_log1p"
    changes <- c(changes, "Applied log1p due to skew > 1; post-hoc back-transformed.") }
    else if (!is.na(shapiro_p) && shapiro_p < alpha) {
      changes <- c(changes, "Normality flagged but log1p not applied (failed refit or negative values).")
    }
  } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
    changes <- c(changes, "Normality flagged; retained scale (no positive-only transform justified).")
  }
  
  # Optional heteroskedastic modeling if Levene < alpha
  selected <- base_fit
  if (!is.na(levene_p) && levene_p < alpha) {
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * time + time:E_phase_trt"))
    fit_try <- function(w) tryCatch(nlme::lme(fixed=form_fixed, random=~1|bottle, data=d0,
                                              method="REML", weights=w,
                                              control=nlme::lmeControl(msMaxIter=200)),
                                    error=function(e) NULL)
    cands <- list(
      none = fit_try(NULL),
      byT  = fit_try(nlme::varIdent(~1|time)),
      byE  = fit_try(nlme::varIdent(~1|E_phase)),
      byS  = fit_try(nlme::varIdent(~1|species)),
      byA  = fit_try(nlme::varIdent(~1|sal_adapt))
    )
    ok <- cands[!vapply(cands, is.null, logical(1))]
    if (length(ok)) {
      aics <- sapply(ok, AIC); pick <- names(which.min(aics))
      selected <- ok[[pick]]
      changes <- c(changes, paste0("Modeled heteroskedasticity via nlme::lme + varIdent (", pick, ") by AIC."))
    }
  }
  
  # Inference engine (drop RE if singular)
  used_engine <- "lmer"
  fit_inf <- tryCatch(refit_lmer_phys(d0, y_use), error=function(e) NULL)
  diag <- quick_diag(fit_inf, data = d0,
                     label = paste0(label, " (DMS(P)t @ t0&t24; E_phase treatment coding)"))
  lines <- c(lines, "Diagnostics (auto, no plots):", diag$lines)
  
  if (is.null(fit_inf)) return(c(lines, "Skipped: inference refit failed.\n"))
  if (inherits(fit_inf, "merMod") && lme4::isSingular(fit_inf, tol=1e-4)) {
    fit_inf <- stats::lm(as.formula(paste(y_use, "~ species * sal_adapt * time + time:E_phase_trt")), data = d0)
    used_engine <- "lm"
    changes <- c(changes, "Random intercept singular; dropped RE and used OLS (lm) for inference.")
  }
  
  aov_tab <- suppressWarnings(car::Anova(fit_inf, type="III"))
  aov_df  <- as.data.frame(aov_tab)
  pcol <- intersect(colnames(aov_df), c("Pr(>F)","Pr(>Chisq)","Pr(>Chi)","Pr..F."))
  if (!length(pcol)) return(c(lines, "Type-III ANOVA failed (no p-value column).", ""))
  
  pullp <- function(term){
    row <- aov_df[rownames(aov_df)==term, , drop=FALSE]
    if (!nrow(row)) return(NA_real_)
    as.numeric(row[[pcol[1]]][1])
  }
  pS  <- pullp("species")
  pA  <- pullp("sal_adapt")
  pT  <- pullp("time")
  pSA <- pullp("species:sal_adapt")
  pST <- pullp("species:time")
  pAT <- pullp("sal_adapt:time")
  pTE <- pullp("time:E_phase_trt")
  
  R2 <- suppressWarnings(tryCatch(MuMIn::r.squaredGLMM(selected), error=function(e) c(NA_real_, NA_real_)))
  
  lines <- c(lines,
             "Assumption checks:", assumptions,
             paste0("Assumptions flagged (α=", alpha, "): ",
                    if (length(flagged)) paste(flagged, collapse = "; ") else "none."),
             if (length(changes)) c("Changes made:", paste0("  - ", changes)) else "Changes made: none.",
             paste0("Inference engine: ", if (used_engine=="lm") "OLS (lm) — random intercept dropped"
                    else "lmer (random-intercept)"),
             "Type-III ANOVA (on inference engine above):",
             paste0("  - species: p=", fmt_p(pS),
                    "; sal_adapt: p=", fmt_p(pA),
                    "; time: p=", fmt_p(pT),
                    "; time×E_phase: p=", fmt_p(pTE)),
             paste0("  - two-way (no aliasing): S×A p=", fmt_p(pSA),
                    "; S×T p=", fmt_p(pST),
                    "; A×T p=", fmt_p(pAT)),
             paste0("Model R² (selected fit) — marginal=", fmt_n(as.numeric(R2[1])),
                    ", conditional=", fmt_n(as.numeric(R2[2]))))
  
  # ----- Gated post-hoc -----
  proceed <- any(c(pS,pA,pT,pSA,pST,pAT,pTE) < alpha, na.rm = TRUE)
  if (!proceed) return(c(lines, "Post-hoc: skipped (no omnibus terms significant at α=0.05).", ""))
  
  trn <- if (!is.null(transform) && transform=="log1p") "log1p" else NULL
  regrid_if <- function(em) if (!is.null(trn) && trn=="log1p") emmeans::regrid(em, transform="response") else em
  as_lines <- function(df){
    if (!nrow(df)) return(character())
    sprintf("  - %s: %s [%s, %s]; p=%s",
            df$contrast %||% df$comparison %||% df$levels,
            fmt_n(df$estimate), fmt_n(df$lower.CL), fmt_n(df$upper.CL), fmt_p(df$p.value))
  }
  
  posthoc_lines <- character()
  
  # Time within species × acclimation (collapse over exposure)
  if (!is.na(pT) && pT < alpha || !is.na(pST) && pST < alpha || !is.na(pAT) && pAT < alpha) {
    emT <- regrid_if(emmeans::emmeans(fit_inf, ~ time, by = c("species","sal_adapt")))
    cmpT <- summary(emmeans::contrast(emT, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Time within species × acclimation):", as_lines(as.data.frame(cmpT)))
  }
  
  # Species within acclimation × time
  if (!is.na(pS) && pS < alpha || !is.na(pST) && pST < alpha) {
    emS <- regrid_if(emmeans::emmeans(fit_inf, ~ species, by = c("sal_adapt","time")))
    cmpS <- summary(emmeans::contrast(emS, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Species within acclimation × time):", as_lines(as.data.frame(cmpS)))
  }
  
  # Acclimation within species × time
  if (!is.na(pA) && pA < alpha || !is.na(pAT) && pAT < alpha) {
    emA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt, by = c("species","time")))
    cmpA <- summary(emmeans::contrast(emA, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Acclimation within species × time):", as_lines(as.data.frame(cmpA)))
  }
  
  # Exposure at t24 only (E25 vs E45) within species × acclimation
  if (!is.na(pTE) && pTE < alpha) {
    emE <- regrid_if(emmeans::emmeans(fit_inf, ~ E_phase_trt, by = c("species","sal_adapt"),
                                      at = list(time = "t24")))
    cmpE <- summary(emmeans::contrast(emE, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    # keep only E25 vs E45 (NoExp is reference & absent at t24)
    cmpE_df <- as.data.frame(cmpE) %>% dplyr::filter(grepl("E25 - E45|E45 - E25", contrast))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Exposure at t24; within species × acclimation):",
                       as_lines(cmpE_df))
  }
  
  c(lines, if (length(posthoc_lines)) posthoc_lines else "Post-hoc: none produced.", "")
}

# -------- Run physiology block and append to the same TXT --------
phys_prose <- c("",
                "==== PHYSIOLOGY (t0 vs t24; exposure only at t24) ====",
                "Physiology metrics are analyzed on raw columns with a time factor (t0 & t24).",
                "At t0 there was no exposure; E_phase uses treatment coding with NoExp as reference.",
                "")

for (y in phys_metrics) {
  if (!y %in% names(d_phys)) {
    phys_prose <- c(phys_prose, hdr(phys_labels[[y]]), "Skipped: column not present.\n"); next
  }
  blk <- tryCatch(analyze_phys_one(d_phys, y, phys_labels[[y]]),
                  error=function(e) c(hdr(phys_labels[[y]]), paste("Analysis failed:", conditionMessage(e)), ""))
  phys_prose <- c(phys_prose, blk)
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path("results", "ANOVA_prose_results.txt")
cat(phys_prose, file = out_txt, sep = "\n", append = TRUE)
message("Physiology results appended for: ", paste(phys_metrics, collapse = ", "))

# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
  library(lme4); library(nlme); library(emmeans); library(car); library(MuMIn)
  library(broom); library(e1071)
})

# ---------- options, seed & guards ----------
set.seed(1)  # reproducibility of any stochastic optimizer steps
.old_contr <- options("contrasts")
options(contrasts = c("contr.sum","contr.poly"))
on.exit(do.call(options, .old_contr), add = TRUE)

if (exists("random", envir = .GlobalEnv) && !is.function(get("random", envir = .GlobalEnv))) {
  try(rm("random", envir = .GlobalEnv), silent = TRUE)
}

# Toggle: if FALSE, avoid keeping singular lmer fits (will try lme/gls fallback for model selection)
ALLOW_SINGULAR <- TRUE
alpha <- 0.05

# ---------- helpers ----------
`%||%` <- function(a,b) if(!is.null(a)) a else b
fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 1e-3, formatC(p, format="e", digits=2),
                                   sprintf("%.3f", p)))
fmt_n <- function(x, d=3) ifelse(is.na(x), "NA", formatC(x, format="f", digits=d))
hdr <- function(txt){ c("", paste0("## ", txt), strrep("-", nchar(txt)+3)) }

# ---------- load & harmonise ----------
chr <- read_csv("clean_chr.csv", show_col_types = FALSE)
pol <- tryCatch(read_csv("clean_pol_modified.csv", show_col_types = FALSE),
                error = function(e) read_csv("clean_pol.csv", show_col_types = FALSE))

dat_all <- bind_rows(chr, pol) %>%
  mutate(
    time_num = suppressWarnings(as.integer(as.character(time))),
    time_std = if_else(time_num == 0L, "t0", "t24"),
    time     = factor(time_std, levels = c("t0","t24")),
    species   = factor(species),
    sal_adapt = factor(sal_adapt, levels = c(25,45), labels = c("A25","A45")),
    sal_exp   = factor(sal_exp,   levels = c(25,45), labels = c("E25","E45")),
    repl      = factor(repl),
    bottle    = interaction(species, sal_adapt, repl, drop = TRUE)
  ) %>%
  select(-time_num)

# ---------- BUILD normalized_df (deterministic t24 sourcing) ----------
keys <- c("species","sal_adapt","sal_exp","repl")

# 1) DMSPp-only metrics (t24) — tolerate optional columns
p24 <- dat_all %>%
  filter(time_std == "t24", id == "DMSPp") %>%
  select(
    all_of(keys),
    mu_POC = u_poc_h_1,               # μ-POC
    u_dmsp_u_poc,                     # μ-DMSP / μ-POC
    dmsp_c_poc_mol_mol,               # DMSP-C : POC (mol:mol)
    any_of(c(                         # OPTIONAL fields
      "up_dmspp",                    # uptake % of DMSPp
      "dmsp_uptake_of_total_dmsp",   # fraction/%
      "dms_puptake_c_poc_mol_mol"    # uptake-C : POC (mol:mol)
    ))
  )

# 2) DMS(P)t-only metrics (t24): direct uptake & loss
t24 <- dat_all %>%
  filter(time_std == "t24", id == "DMS(P)t") %>%
  select(
    all_of(keys),
    uptake_alt = d3_p_taken_up,
    loss_alt   = d3_p_lost_demethylated
  )

# 3) μ-DMSPd / μ-POC: numerator from DMSPd (t24), denominator μ-POC from DMSPp (t24)
d24 <- dat_all %>%
  filter(time_std == "t24", id == "DMSPd") %>%
  select(all_of(keys), mu_DMSPd = u_dmsp_h_1)

mu_ratio_tbl <- p24 %>%
  select(all_of(keys), mu_POC) %>%
  left_join(d24, by = keys) %>%
  mutate(muDMSPd_over_muPOC = if_else(
    is.finite(mu_DMSPd) & is.finite(mu_POC) & mu_POC != 0,
    mu_DMSPd / mu_POC, NA_real_
  )) %>%
  select(all_of(keys), muDMSPd_over_muPOC)

# 4) Assemble normalized_df (all rows are t24 + tagged)
normalized_df <- p24 %>%
  left_join(t24,          by = keys) %>%
  left_join(mu_ratio_tbl, by = keys) %>%
  mutate(
    time_std = "t24",
    time     = factor("t24", levels = c("t0","t24")),
    id       = "Normalized"
  ) %>%
  select(
    all_of(keys), time, time_std, id,
    mu_POC, u_dmsp_u_poc, dmsp_c_poc_mol_mol,
    any_of(c("up_dmspp", "dmsp_uptake_of_total_dmsp", "dms_puptake_c_poc_mol_mol")),
    muDMSPd_over_muPOC,
    uptake_alt, loss_alt
  )

# ---------- analyses (scoped to t24 inside analyze_one) ----------
analyses <- list(
  # Production (raw, by source id)
  list(name="mu_POC_DMSPp",           dat=dat_all,       y="u_poc_h_1",              id="DMSPp",   label="μ-POC (from DMSPp)"),
  list(name="mu_DMSP_from_DMSPp",     dat=dat_all,       y="u_dmsp_h_1",             id="DMSPp",   label="μ-DMSP (from DMSPp)"),
  list(name="mu_DMSP_from_DMS",       dat=dat_all,       y="u_dmsp_h_1",             id="DMS",     label="μ-DMSP (from DMS)"),
  list(name="uptake_d3",              dat=dat_all,       y="d3_p_taken_up",          id="DMS(P)t", label="D3-DMSP taken up"),
  list(name="loss_d3",                dat=dat_all,       y="d3_p_lost_demethylated", id="DMS(P)t", label="D3-DMSP lost (demethylated)"),
  
  # Normalized (t24 bundle, already t24-scoped by construction)
  list(name="norm_u_dmsp_u_poc",      dat=normalized_df, y="u_dmsp_u_poc",           id="Normalized", label="μ-DMSP / μ-POC (direct)"),
  list(name="norm_uptake_alt_direct", dat=normalized_df, y="uptake_alt",             id="Normalized", label="Uptake ALT (direct d3_p_taken_up)"),
  list(name="norm_muDMSPd_muPOC",     dat=normalized_df, y="muDMSPd_over_muPOC",     id="Normalized", label="μ-DMSPd / μ-POC"),
  list(name="norm_dmsp_c_poc",        dat=normalized_df, y="dmsp_c_poc_mol_mol",     id="Normalized", label="DMSP-C : POC (mol:mol)"),
  
  # EXACT SCOPES (t24 + correct source IDs). Will auto-skip if column missing.
  list(name="up_dmspp_t24",                  dat=dat_all, y="up_dmspp",                  id="DMSPp",   label="Uptake (% of DMSPp)"),
  list(name="dmsp_uptake_total_t24",         dat=dat_all, y="dmsp_uptake_of_total_dmsp", id="DMSPp",   label="DMSP uptake of total DMSP"),
  list(name="dms_puptake_c_poc_t24",         dat=dat_all, y="dms_puptake_c_poc_mol_mol", id="DMSPp",   label="DMSP uptake-C : POC (mol:mol)")
)

# ---------- core engines ----------
refit_lmer <- function(data, resp){
  form <- as.formula(paste(resp, "~ species * sal_adapt * sal_exp + (1|bottle)"))
  lme4::lmer(form, data = data, REML = TRUE,
             control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

analyze_one <- function(dat, y, label, id_filter, time_filter = "t24"){
  lines <- c(hdr(label))
  
  # filter scope
  d0 <- dat %>%
    filter(time_std == time_filter, id == id_filter) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp, repl) %>%
    mutate(bottle = interaction(species, sal_adapt, repl, drop = TRUE)) %>%
    mutate(.y = suppressWarnings(as.numeric(.data[[y]]))) %>%
    filter(is.finite(.y))
  
  # structural guards
  if (nrow(d0) < 8L ||
      n_distinct(d0$species)   < 2 ||
      n_distinct(d0$sal_adapt) < 2 ||
      n_distinct(d0$sal_exp)   < 2) {
    return(c(lines, "Skipped: insufficient data or factors lack ≥2 levels.\n"))
  }
  if (length(unique(d0$.y)) < 2) {
    return(c(lines, "Skipped: response has < 2 unique finite values.\n"))
  }
  
  # Baseline lmer & diagnostics
  base_fit <- tryCatch(refit_lmer(d0, ".y"), error=function(e) NULL)
  if (is.null(base_fit)) return(c(lines, "Skipped: baseline lmer failed.\n"))
  
  r <- residuals(base_fit, type="pearson")
  shapiro_p <- tryCatch(shapiro.test(as.numeric(r))$p.value, error=function(e) NA_real_)
  grp <- with(d0, interaction(species, sal_adapt, sal_exp, drop = TRUE))
  levene_p <- tryCatch(car::leveneTest(r ~ grp, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
  
  assumptions <- c(
    paste0("Normality (Shapiro–Wilk on residuals): p = ", fmt_p(shapiro_p)),
    paste0("Homoskedasticity (Levene across S×A×E): p = ", fmt_p(levene_p))
  )
  flagged <- c()
  if (!is.na(shapiro_p) && shapiro_p < alpha) flagged <- c(flagged, paste0("Normality (p=", fmt_p(shapiro_p), ")"))
  if (!is.na(levene_p)  && levene_p  < alpha) flagged <- c(flagged,  paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
  
  changes <- character()
  transform <- NULL
  y_use <- ".y"
  
  # Optional transform for strong positive skew
  sk <- suppressWarnings(e1071::skewness(d0$.y, na.rm = TRUE))
  if (is.finite(sk) && sk > 1 && all(d0$.y >= 0, na.rm = TRUE)) {
    d0$.y_log1p <- log1p(d0$.y)
    base_fit2 <- tryCatch(refit_lmer(d0, ".y_log1p"), error=function(e) NULL)
    if (!is.null(base_fit2)) {
      base_fit <- base_fit2; transform <- "log1p"; y_use <- ".y_log1p"
      changes <- c(changes, "Applied log1p due to skew > 1; post-hoc back-transformed.")
      r2 <- residuals(base_fit, type="pearson")
      levene_p <- tryCatch(car::leveneTest(r2 ~ grp, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
      flagged <- flagged[!grepl("^Homoskedasticity", flagged)]
      if (!is.na(levene_p) && levene_p < alpha)
        flagged <- c(flagged, paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
    } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
      changes <- c(changes, "Normality flagged but log1p not applied (failed refit or negative values).")
    }
  } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
    changes <- c(changes, "Normality flagged; retained scale (no positive-only transform justified).")
  }
  
  # Heteroskedastic candidates via nlme::lme with varIdent if Levene < alpha (AIC selection)
  selected <- base_fit
  if (!is.na(levene_p) && levene_p < alpha) {
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * sal_exp"))
    fit_try <- function(w) tryCatch(nlme::lme(fixed=form_fixed, random=~1|bottle, data=d0,
                                              method="REML", weights=w, control=nlme::lmeControl(msMaxIter=200)),
                                    error=function(e) NULL)
    cands <- list(
      none = fit_try(NULL),
      byE  = fit_try(nlme::varIdent(~1|sal_exp)),
      byS  = fit_try(nlme::varIdent(~1|species)),
      byA  = fit_try(nlme::varIdent(~1|sal_adapt)),
      byAE = fit_try(nlme::varIdent(~1|sal_adapt:sal_exp))
    )
    ok <- cands[!vapply(cands, is.null, logical(1))]
    if (length(ok)) {
      aics <- sapply(ok, AIC)
      pick <- names(which.min(aics))
      selected <- ok[[pick]]
      changes <- c(changes, paste0("Modeled heteroskedasticity via nlme::lme + varIdent (", pick, ") by AIC."))
    }
  }
  
  # Optionally avoid singular lmer for selection block
  if (!ALLOW_SINGULAR && inherits(selected, "merMod") && lme4::isSingular(selected, tol=1e-4)) {
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * sal_exp"))
    sel2 <- tryCatch(nlme::lme(fixed=form_fixed, random=~1|bottle, data=d0, method="REML",
                               control=nlme::lmeControl(msMaxIter=200)), error=function(e) NULL)
    if (is.null(sel2)) sel2 <- tryCatch(nlme::gls(form_fixed, data=d0, method="REML"), error=function(e) NULL)
    if (!is.null(sel2)) { selected <- sel2; changes <- c(changes, "Avoided singular lmer by switching to lme/gls for selection.") }
  }
  
  # ---------- INFERENCE ENGINE (lmer unless singular → lm) ----------
  used_engine <- "lmer"
  rand_var <- NA_real_
  if (inherits(base_fit, "merMod")) {
    vc_tab <- tryCatch(as.data.frame(VarCorr(base_fit)), error=function(e) NULL)
    if (!is.null(vc_tab)) {
      rv <- tryCatch(vc_tab %>% dplyr::filter(grp=="bottle", var1=="(Intercept)") %>% dplyr::pull(vcov), error=function(e) numeric())
      rand_var <- rv[1] %||% NA_real_
    }
  }
  
  fit_inf <- tryCatch(refit_lmer(d0, y_use), error=function(e) NULL)
  if (is.null(fit_inf)) return(c(lines, "Skipped: inference refit failed.\n"))
  if (inherits(fit_inf, "merMod") && lme4::isSingular(fit_inf, tol=1e-4)) {
    # Drop random intercept and use OLS for inference on fixed effects
    fit_inf <- stats::lm(as.formula(paste(y_use, "~ species * sal_adapt * sal_exp")), data = d0)
    used_engine <- "lm"
    changes <- c(changes,
                 sprintf("Random intercept singular (Var[bottle] ≈ %s); dropped RE and used OLS (lm) for inference.",
                         fmt_n(rand_var, d=6)))
  }
  
  # Type-III ANOVA on inference engine
  aov_tab <- suppressWarnings(car::Anova(fit_inf, type="III"))
  aov_df  <- as.data.frame(aov_tab)
  pcol <- intersect(colnames(aov_df), c("Pr(>F)","Pr(>Chisq)","Pr(>Chi)","Pr..F."))
  if (!length(pcol)) return(c(lines, "Type-III ANOVA failed (no p-value column).", ""))
  
  getp <- function(term){
    row <- aov_df[rownames(aov_df)==term, , drop=FALSE]
    if (!nrow(row)) return(NA_real_)
    as.numeric(row[[pcol[1]]][1])
  }
  p3  <- getp("species:sal_adapt:sal_exp")
  pSA <- getp("species:sal_adapt")
  pSE <- getp("species:sal_exp")
  pAE <- getp("sal_adapt:sal_exp")
  pS  <- getp("species")
  pA  <- getp("sal_adapt")
  pE  <- getp("sal_exp")
  
  # R2 from selected model (AIC-chosen; may be lmer/lme/gls). This is for description only.
  R2 <- suppressWarnings(tryCatch(MuMIn::r.squaredGLMM(selected), error=function(e) c(NA_real_, NA_real_)))
  
  # --- Prose reporting up to ANOVA ---
  lines <- c(lines,
             "Assumption checks:",
             paste0("  - ", assumptions),
             paste0("Assumptions flagged (α=", alpha, "): ",
                    if (length(flagged)) paste(flagged, collapse = "; ") else "none."),
             if (length(changes)) c("Changes made:", paste0("  - ", changes)) else "Changes made: none.",
             paste0("Inference engine: ",
                    if (used_engine=="lm") "OLS (lm) — random intercept dropped due to singularity"
                    else "lmer (random-intercept)"),
             if (!is.na(rand_var)) paste0("  - Estimated Var[bottle] from lmer: ", fmt_n(rand_var, d=6)) else NULL,
             "Type-III ANOVA (on inference engine above):",
             paste0("  - species: p=", fmt_p(pS),
                    "; sal_adapt: p=", fmt_p(pA),
                    "; sal_exp: p=", fmt_p(pE)),
             paste0("  - species:sal_adapt: p=", fmt_p(pSA),
                    "; species:sal_exp: p=", fmt_p(pSE),
                    "; sal_adapt:sal_exp: p=", fmt_p(pAE)),
             paste0("  - species:sal_adapt:sal_exp: p=", fmt_p(p3)),
             paste0("Model R² (selected fit) — marginal=", fmt_n(as.numeric(R2[1])),
                    ", conditional=", fmt_n(as.numeric(R2[2]))))
  
  # --- GATED post-hoc: proceed only if omnibus indicates it ---
  proceed <- any(c(p3,pSA,pSE,pAE,pS,pA,pE) < alpha, na.rm = TRUE)
  if (!proceed) {
    lines <- c(lines, "Post-hoc: skipped (no omnibus terms significant at α=0.05).", "")
    return(lines)
  }
  
  # Choose what to compare:
  do_all <- any(c(p3,pSA,pSE,pAE) < alpha, na.rm = TRUE)
  trn <- if (!is.null(transform) && transform=="log1p") "log1p" else NULL
  
  posthoc_lines <- character()
  as_lines <- function(df){
    if (!nrow(df)) return(character())
    sprintf("  - %s: %s [%s, %s]; p=%s",
            df$contrast %||% df$comparison %||% df$levels,
            fmt_n(df$estimate), fmt_n(df$lower.CL), fmt_n(df$upper.CL), fmt_p(df$p.value))
  }
  regrid_if <- function(em) if (!is.null(trn) && trn=="log1p") emmeans::regrid(em, transform="response") else em
  
  if (do_all) {
    em1 <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_exp, by = c("species","sal_adapt")))
    cmp1 <- summary(emmeans::contrast(em1, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (E within species × acclimation):", as_lines(as.data.frame(cmp1)))
    
    emS <- regrid_if(emmeans::emmeans(fit_inf, ~ species, by = c("sal_adapt","sal_exp")))
    cmpS <- summary(emmeans::contrast(emS, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Species within acclimation × exposure):", as_lines(as.data.frame(cmpS)))
    
    emA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt, by = c("species","sal_exp")))
    cmpA <- summary(emmeans::contrast(emA, "pairwise", adjust="holm"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (Acclimation within species × exposure):", as_lines(as.data.frame(cmpA)))
    
    em_all <- regrid_if(emmeans::emmeans(fit_inf, ~ species * sal_adapt * sal_exp))
    cmp_all <- summary(pairs(em_all, adjust="tukey"), infer=c(TRUE,TRUE))
    posthoc_lines <- c(posthoc_lines, "Post-hoc (All 3-way cell means; Tukey):", as_lines(as.data.frame(cmp_all)))
  } else {
    if (!is.na(pS) && pS < alpha) {
      em_mS <- regrid_if(emmeans::emmeans(fit_inf, ~ species))
      cmp_mS <- summary(pairs(em_mS, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Species; marginal):", as_lines(as.data.frame(cmp_mS)))
    }
    if (!is.na(pA) && pA < alpha) {
      em_mA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt))
      cmp_mA <- summary(pairs(em_mA, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Acclimation; marginal):", as_lines(as.data.frame(cmp_mA)))
    }
    if (!is.na(pE) && pE < alpha) {
      em_mE <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_exp))
      cmp_mE <- summary(pairs(em_mE, adjust="holm"), infer=c(TRUE,TRUE))
      posthoc_lines <- c(posthoc_lines, "Post-hoc (Exposure; marginal):", as_lines(as.data.frame(cmp_mE)))
    }
  }
  
  c(lines, if (length(posthoc_lines)) posthoc_lines else "Post-hoc: none produced.", "")
}

# ---------- Run all analyses & write TXT ----------
prose <- c("==== MIXED-MODEL RESULTS (prose) ====",
           paste0("Alpha = ", alpha),
           "Type-III ANOVA computed on a clean lmer refit; if the random intercept is singular,",
           "we drop it and use OLS (lm) for inference on fixed effects. Post-hoc uses the same engine.",
           "If Levene < alpha, model selection considers nlme::lme + varIdent by AIC (emmeans not run on lme).",
           "")

for (a in analyses) {
  if (!a$y %in% names(a$dat)) {
    prose <- c(prose, hdr(a$label), "Skipped: response column not present.\n")
    next
  }
  blk <- tryCatch(analyze_one(a$dat, a$y, a$label, a$id, time_filter = "t24"),
                  error=function(e) c(hdr(a$label), paste("Analysis failed:", conditionMessage(e)), ""))
  prose <- c(prose, blk)
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path("results", "ANOVA_prose_results.txt")
writeLines(prose, out_txt)
cat("Wrote TXT to: ", normalizePath(out_txt, winslash = "/"), "\n", sep="")

# ---------- (LaTeX export intentionally disabled) ----------
# library(xtable)
# make_family_table_tex_all <- function(...) { }
# combined_base <- file.path("results", "ALL_contrasts_combined")
# make_family_table_tex_all(ok_names, results, combined_base)

