# ============================================================
# ONE-BLOCK PIPELINE (prose TXT only; no LaTeX)
# t24-only changes: use raw columns directly by ID (no re-compute of responses)
# ============================================================
set.seed(42) #this is the first comment

setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")

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
# PHYSIOLOGY — t0→t24 DELTAS (paired by bottle; exposure only at t24)
# Consistent with t24-only blocks: model Δ ~ species * sal_adapt * sal_exp + (1|bottle)
# Baselines (t0) from id == "DMS(P)t"; exposure applied only at t24 branches.
# Deltas computed for: mean cell size, cell density, Fv/Fm ("Yield").
# ============================================================

# Build paired dataset per bottle branch (control/shock at t24)
d_phys_src <- dat_all %>%
  dplyr::filter(id == "DMS(P)t") %>%
  mutate(bottle = interaction(species, sal_adapt, repl, drop = TRUE))

# Baseline at t0 (optionally scaled density to bottle start, as in earlier runs)
base_t0 <- d_phys_src %>%
  dplyr::filter(time == "t0") %>%
  mutate(dens_t0_bottle = 0.6 * cell_dens_avr) %>%   # keep prior convention
  dplyr::select(
    species, sal_adapt, repl,
    yield_t0 = yield,
    dens_t0  = dens_t0_bottle,
    size_t0  = cell_size_avr
  )

# t24 observed (exposure present here)
t24 <- d_phys_src %>%
  dplyr::filter(time == "t24") %>%
  dplyr::select(
    species, sal_adapt, sal_exp, repl,
    yield_t24 = yield,
    dens_t24  = cell_dens_avr,
    size_t24  = cell_size_avr
  )

# Pair & compute deltas (and simple % change for quick summaries if needed)
phys_delta <- t24 %>%
  dplyr::inner_join(base_t0, by = c("species","sal_adapt","repl")) %>%
  dplyr::mutate(
    bottle   = interaction(species, sal_adapt, repl, drop = TRUE),
    d_yield  = yield_t24 - yield_t0,
    d_dens   = dens_t24  - dens_t0,
    d_size   = size_t24  - size_t0,
    pct_yield = ifelse(is.finite(yield_t0) & yield_t0 != 0,
                       (yield_t24 / yield_t0 - 1) * 100, NA_real_),
    pct_dens  = ifelse(is.finite(dens_t0)  & dens_t0  != 0,
                       (dens_t24  / dens_t0  - 1) * 100, NA_real_),
    pct_size  = ifelse(is.finite(size_t0)  & size_t0  != 0,
                       (size_t24  / size_t0  - 1) * 100, NA_real_)
  )

# Helper: analyze one Δ-metric with the SAME output structure as other endpoints
analyze_phys_delta <- function(dat, y, nice_label){
  lines <- c(hdr(nice_label),
             "Paired t0→t24 analysis on deltas (Δ = t24 − t0).",
             "Model: Δ ~ species × sal_adapt × sal_exp + (1|bottle); exposure only at t24.")
  
  d0 <- dat %>%
    tidyr::drop_na(species, sal_adapt, sal_exp, repl, bottle) %>%
    dplyr::mutate(.y = suppressWarnings(as.numeric(.data[[y]]))) %>%
    dplyr::filter(is.finite(.y))
  
  if (nrow(d0) < 8L ||
      dplyr::n_distinct(d0$species)   < 2 ||
      dplyr::n_distinct(d0$sal_adapt) < 2 ||
      dplyr::n_distinct(d0$sal_exp)   < 2) {
    return(c(lines, "Skipped: insufficient data or factors lack ≥2 levels.\n"))
  }
  if (length(unique(d0$.y)) < 2) return(c(lines, "Skipped: response has < 2 unique finite values.\n"))
  
  # Base fit (random-intercept)
  base_fit <- tryCatch(refit_lmer(d0, ".y"), error=function(e) NULL)
  if (is.null(base_fit)) return(c(lines, "Skipped: baseline lmer failed.\n"))
  
  # Diagnostics (grouping by S×A×E to mirror other t24 blocks)
  r <- residuals(base_fit, type="pearson")
  shapiro_p <- tryCatch(shapiro.test(as.numeric(r))$p.value, error=function(e) NA_real_)
  grp <- with(d0, interaction(species, sal_adapt, sal_exp, drop = TRUE))
  levene_p <- tryCatch(car::leveneTest(r ~ grp, center = median)[[1,"Pr(>F)"]], error=function(e) NA_real_)
  
  assumptions <- c(
    paste0("  - Normality (Shapiro–Wilk on residuals): p = ", fmt_p(shapiro_p)),
    paste0("  - Homoskedasticity (Levene across S×A×E): p = ", fmt_p(levene_p))
  )
  flagged <- c()
  if (!is.na(shapiro_p) && shapiro_p < alpha) flagged <- c(flagged, paste0("Normality (p=", fmt_p(shapiro_p), ")"))
  if (!is.na(levene_p)  && levene_p  < alpha) flagged <- c(flagged,  paste0("Homoskedasticity (p=", fmt_p(levene_p), ")"))
  
  # Transform only if positive & strongly skewed (deltas can be negative → usually no transform)
  changes <- character(); y_use <- ".y"; transform <- NULL
  sk <- suppressWarnings(e1071::skewness(d0$.y, na.rm = TRUE))
  if (is.finite(sk) && sk > 1 && all(d0$.y >= 0, na.rm = TRUE)) {
    d0$.y_log1p <- log1p(d0$.y)
    base_fit2 <- tryCatch(refit_lmer(d0, ".y_log1p"), error=function(e) NULL)
    if (!is.null(base_fit2)) { base_fit <- base_fit2; y_use <- ".y_log1p"; transform <- "log1p"
    changes <- c(changes, "Applied log1p due to skew > 1; post-hoc back-transformed.") }
  } else if (!is.na(shapiro_p) && shapiro_p < alpha) {
    changes <- c(changes, "Normality flagged; retained scale (Δ allows negatives).")
  }
  
  # Optional heteroskedastic modeling (for R² reporting)
  selected <- base_fit
  if (!is.na(levene_p) && levene_p < alpha) {
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * sal_exp"))
    fit_try <- function(w) tryCatch(nlme::lme(fixed=form_fixed, random=~1|bottle, data=d0,
                                              method="REML", weights=w,
                                              control=nlme::lmeControl(msMaxIter=200)),
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
  
  # Inference engine (drop RE if singular) — identical to other endpoints
  fit_inf <- tryCatch(refit_lmer(d0, y_use), error=function(e) NULL)
  lines <- c(lines, "Diagnostics (auto, no plots):",
             quick_diag(fit_inf, data = d0, label = paste0(nice_label, " (Δ @ t24; paired to t0)"))$lines)
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
             paste0("  - Normality (Shapiro–Wilk on residuals): p = ", fmt_p(shapiro_p)),
             paste0("  - Homoskedasticity (Levene across S×A×E): p = ", fmt_p(levene_p)),
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
  
  # Post-hoc (same families as other t24 endpoints)
  regrid_if <- function(em) if (!is.null(transform) && transform=="log1p") emmeans::regrid(em, transform="response") else em
  as_lines <- function(df) if (!nrow(df)) character() else
    sprintf("  - %s: %s [%s, %s]; p=%s",
            df$contrast, fmt_n(df$estimate), fmt_n(df$lower.CL), fmt_n(df$upper.CL), fmt_p(df$p.value))
  
  posthoc_lines <- character()
  emE <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_exp,   by = c("species","sal_adapt")))
  emS <- regrid_if(emmeans::emmeans(fit_inf, ~ species,   by = c("sal_adapt","sal_exp")))
  emA <- regrid_if(emmeans::emmeans(fit_inf, ~ sal_adapt, by = c("species","sal_exp")))
  em3 <- regrid_if(emmeans::emmeans(fit_inf, ~ species * sal_adapt * sal_exp))
  
  cmpE   <- summary(emmeans::contrast(emE, "pairwise", adjust="holm"),  infer=c(TRUE,TRUE))
  cmpS   <- summary(emmeans::contrast(emS, "pairwise", adjust="holm"),  infer=c(TRUE,TRUE))
  cmpA   <- summary(emmeans::contrast(emA, "pairwise", adjust="holm"),  infer=c(TRUE,TRUE))
  cmpAll <- summary(pairs(em3, adjust="tukey"),                          infer=c(TRUE,TRUE))
  
  posthoc_lines <- c(posthoc_lines,
                     "Post-hoc (E within species × acclimation):",          as_lines(as.data.frame(cmpE)),
                     "Post-hoc (Species within acclimation × exposure):",   as_lines(as.data.frame(cmpS)),
                     "Post-hoc (Acclimation within species × exposure):",   as_lines(as.data.frame(cmpA)),
                     "Post-hoc (All 3-way cell means; Tukey):",             as_lines(as.data.frame(cmpAll)))
  c(lines, posthoc_lines, "")
}

# Which Δ metrics to analyse
phys_delta_metrics <- c(
  d_size  = "Physiology Δ: Mean cell size",
  d_dens  = "Physiology Δ: Cell density",
  d_yield = "Physiology Δ: Yield (Fv/Fm)"
)

# Run and append to the SAME TXT as before
phys_delta_prose <- c("",
                      "==== PHYSIOLOGY (t0→t24 deltas; exposure only at t24) ====",
                      "Each t24 branch (E25/E45) is paired to its own bottle's t0 baseline; deltas are modeled with the same S×A×E fixed-effects structure as other endpoints.",
                      "")

for (nm in names(phys_delta_metrics)) {
  ycol <- nm
  if (!ycol %in% names(phys_delta)) {
    phys_delta_prose <- c(phys_delta_prose, hdr(phys_delta_metrics[[nm]]), "Skipped: column not present.\n")
    next
  }
  blk <- tryCatch(analyze_phys_delta(phys_delta, ycol, phys_delta_metrics[[nm]]),
                  error=function(e) c(hdr(phys_delta_metrics[[nm]]), paste("Analysis failed:", conditionMessage(e)), ""))
  phys_delta_prose <- c(phys_delta_prose, blk)
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path("results", "ANOVA_prose_results.txt")
cat(phys_delta_prose, file = out_txt, sep = "\n", append = TRUE)
message("Physiology Δ results appended for: ", paste(names(phys_delta_metrics), collapse = ", "))
