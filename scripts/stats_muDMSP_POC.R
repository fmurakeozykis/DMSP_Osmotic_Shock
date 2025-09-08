# ---- Packages ----
library(dplyr)
library(readr)
library(car)
library(emmeans)

# ----- add extra mu-POC ----

pol <- read_csv("clean_pol.csv")

# Add replicate C row with only known fields + μ-POC; others NA
new_row <- tibble(
  species    = "Polarella",
  sal_adapt  = 45,
  sal_exp    = 25,
  time       = 24,
  repl       = "C",
  id         = "DMSPp",
  u_poc_h_1  = 0.0029
)

# Bind and save
pol_mod <- bind_rows(pol, new_row)

write_csv(pol_mod, "clean_pol_modified.csv")


# ---- Load & combine ----
chr <- read_csv("clean_chr.csv")
pol <- read_csv("clean_pol_modified.csv")
df  <- bind_rows(chr, pol)

# Factor coding
df <- df %>%
  mutate(
    species   = factor(species),
    sal_adapt = factor(sal_adapt, levels = c(25,45), labels = c("A25","A45")),
    sal_exp   = factor(sal_exp,   levels = c(25,45), labels = c("E25","E45")),
    time      = factor(time, levels = c(0,24), labels = c("t0","t24"))
  )

# ---- Helper function to run ANOVA with assumption checks ----

run_anova <- function(dat, y, label) {
  # 1) Keep rows with non-missing outcome & required factors
  d <- dat %>%
    dplyr::filter(!is.na(.data[[y]])) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp)
  
  # 2) Check replicates per cell
  cell_vars <- c("species","sal_adapt","sal_exp")
  cell_n <- d %>%
    dplyr::count(dplyr::across(dplyr::all_of(cell_vars)), name = "n")
  
  # Decide whether to use means: singletons OR too few total rows
  need_means <- any(cell_n$n < 2) || nrow(d) < 8
  
  if (need_means) {
    # Warn about which cells forced the change (if any)
    bad <- cell_n %>% dplyr::filter(n < 2)
    if (nrow(bad) > 0) {
      msg_cells <- bad %>%
        dplyr::mutate(cell = paste(species, sal_adapt, sal_exp, sep="/")) %>%
        dplyr::pull(cell) %>% paste(collapse = ", ")
      message("Using cell means for ", label, 
              " due to singleton cells (n<2): ", msg_cells)
    } else {
      message("Using cell means for ", label, 
              " due to too few rows (n<8).")
    }
    
    # Collapse to one mean per cell
    d <- d %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cell_vars))) %>%
      dplyr::summarise(!!y := mean(.data[[y]], na.rm = TRUE), .groups = "drop")
  }
  
  # If still too small, exit early
  if (nrow(d) < 8) {
    message("Too few rows for ", label, " even after checking/averaging.")
    return(NULL)
  }
  
  # 3) Fit the 3-way model
  form <- stats::as.formula(paste(y, "~ species * sal_adapt * sal_exp"))
  fit  <- stats::lm(form, data = d)
  
  # 4) Assumption checks
  resids <- residuals(fit)
  if (length(resids) < 3 || stats::sd(resids, na.rm = TRUE) < .Machine$double.eps) {
    shapiro_p <- NA_real_
    cat("Shapiro skipped: residuals have zero/near-zero variance or too few values.\n")
  } else if (length(resids) > 5000) {
    shapiro_p <- NA_real_
  } else {
    shapiro_p <- stats::shapiro.test(resids)$p.value
  }
  
  lev3_p <- tryCatch(
    car::leveneTest(d[[y]] ~ interaction(d$species, d$sal_adapt, d$sal_exp))[[1,"Pr(>F)"]],
    error = function(e) NA_real_
  )
  lev1_p <- tryCatch(
    car::leveneTest(d[[y]] ~ d$sal_exp)[[1,"Pr(>F)"]],
    error = function(e) NA_real_
  )
  ncv_p  <- tryCatch(car::ncvTest(fit)$`p`, error = function(e) NA_real_)
  
  cat("\n============================\n", label, " — Assumption checks\n", sep="")
  cat("Shapiro (residuals) p =", round(shapiro_p, 4), "\n")
  cat("Levene (3-way groups) p =", round(lev3_p, 4), "\n")
  cat("Levene (by sal_exp)  p =", round(lev1_p, 4), "\n")
  cat("NCV test p =", round(ncv_p, 4), "\n")
  
  normal_ok <- is.na(shapiro_p) || shapiro_p > 0.05
  var_ok    <- (is.na(lev3_p) || lev3_p > 0.05) && (is.na(lev1_p) || lev1_p > 0.05)
  ncv_ok    <- is.na(ncv_p)    || ncv_p > 0.05
  use_hc3   <- !(normal_ok && var_ok && ncv_ok)
  
  cat(if (use_hc3) "\nUsing Type-III ANOVA with White HC3 robust SEs\n"
      else          "\nUsing standard Type-III ANOVA\n")
  
  aov_tab <- if (use_hc3) {
    tryCatch(car::Anova(fit, type = 3, white.adjust = "hc3"),
             error = function(e) {
               cat("\nHC3 failed (often leverage=1). Falling back to standard Type-III.\n")
               car::Anova(fit, type = 3)
             })
  } else {
    car::Anova(fit, type = 3)
  }
  print(aov_tab)
  
  # 5) Post-hoc Tukey: E25 vs E45 within species × sal_adapt
  cat("\nPost-hoc: sal_exp within species × sal_adapt\n")
  #emm <- emmeans::emmeans(fit, ~ sal_exp | species + sal_adapt)
  #print(emmeans::pairs(emm, adjust = "tukey"))
  emm <- emmeans(fit3, ~ sal_exp | species * sal_adapt,
                 vcov. = car::hccm(fit3, type = "hc3"))  # keep robust SEs
  pairs(emm, adjust = "tukey") |> summary(infer = c(TRUE, TRUE))
  
  }



# ---- Metrics and IDs ----
metrics_info <- list(
  list(col = "u_poc_h_1", label = "μ-POC (DMSPp)", id_filter = "DMSPp"),
  list(col = "u_dmsp_h_1", label = "μ-DMSP (DMS)", id_filter = "DMS"),
  list(col = "u_dmsp_h_1", label = "μ-DMSP (DMSPp)", id_filter = "DMSPp"),
  list(col = "u_dmsp_h_1", label = "μ-DMSP (DMSPd)", id_filter = "DMSPd")
)

# ---- Run analyses ----
for (info in metrics_info) {
  sub <- df %>% filter(time == "t24", id == info$id_filter)
  if (info$col %in% names(sub)) {
    run_anova2(sub, info$col, info$label)
  }
}



















################


run_anova <- function(dat, y, label) {
  d <- dat %>% filter(!is.na(.data[[y]]),
                      !is.na(species), !is.na(sal_adapt), !is.na(sal_exp))
  if (nrow(d) < 8) {
    message("Too few rows for ", label)
    return(NULL)
  }
  
  form <- as.formula(paste(y, "~ species * sal_adapt * sal_exp"))
  fit  <- lm(form, data = d)
  
  # --- Assumption checks (info only) ---
  resids    <- residuals(fit)
  shapiro_p <- if (length(resids) >= 3 && length(resids) <= 5000) shapiro.test(resids)$p.value else NA_real_
  lev3_p    <- tryCatch(car::leveneTest(d[[y]] ~ interaction(d$species, d$sal_adapt, d$sal_exp))[[1,"Pr(>F)"]],
                        error = function(e) NA_real_)
  lev1_p    <- tryCatch(car::leveneTest(d[[y]] ~ d$sal_exp)[[1,"Pr(>F)"]],
                        error = function(e) NA_real_)
  ncv_p     <- tryCatch(car::ncvTest(fit)$`p`, error = function(e) NA_real_)
  
  cat("\n============================\n", label, " — Assumption checks\n", sep="")
  cat("Shapiro (residuals) p =", round(shapiro_p, 4), "\n")
  cat("Levene (3-way groups) p =", round(lev3_p, 4), "\n")
  cat("Levene (by sal_exp)  p =", round(lev1_p, 4), "\n")
  cat("NCV test p =", round(ncv_p, 4), "\n")
  
  # --- Try robust Type-III (HC3); if it errors (leverage=1), fall back to standard ---
  aov_tab <- tryCatch({
    cat("\nTrying Type-III ANOVA with White HC3 robust SEs ...\n")
    car::Anova(fit, type = 3, white.adjust = "hc3")
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("\nHC3 failed (likely leverage=1 from a single-replicate cell).\n",
        "Falling back to standard Type-III ANOVA. Interpret with caution.\n", sep = "")
    car::Anova(fit, type = 3)
  })
  print(aov_tab)
  
  # --- Post-hoc: sal_exp within species × sal_adapt (still valid for either fit) ---
  cat("\nPost-hoc: sal_exp within species × sal_adapt (Tukey)\n")
  emm <- emmeans::emmeans(fit, ~ sal_exp | species + sal_adapt)
  print(pairs(emm, adjust = "tukey"))
  
  invisible(list(model = fit, anova = aov_tab, emm = emm))
}
# Keep ONLY the second run_anova() you defined above.

# Subset to the right rows
sub_dmspc_poc <- df %>%
  dplyr::filter(time == "t24", id == "DMSPp")

# Run the 3-way ANOVA (species × sal_adapt × sal_exp)
run_anova(sub_dmspc_poc,
          y = "dmsp_c_poc_mol_mol",
          label = "DMSP-C : POC (mol:mol)")

########################

# Force standard Type-III for this run
run_anova_dmspc_poc <- function(dat, y, label) {
  d <- dat %>% filter(!is.na(.data[[y]]),
                      !is.na(species), !is.na(sal_adapt), !is.na(sal_exp))
  
  form <- as.formula(paste(y, "~ species * sal_adapt * sal_exp"))
  fit  <- lm(form, data = d)
  
  cat("\n============================\n", label, "\n", sep="")
  print(car::Anova(fit, type = 3))
  
  cat("\nPost-hoc: sal_exp within species × sal_adapt\n")
  emm <- emmeans::emmeans(fit, ~ sal_exp | species + sal_adapt)
  print(pairs(emm, adjust = "tukey"))
}

# Run
run_anova_dmspc_poc(
  df %>% filter(time == "t24", id == "DMSPp"),
  y = "dmsp_c_poc_mol_mol",
  label = "DMSP-C : POC (mol:mol)"
)




# ---- Dependencies ----
library(dplyr)
library(car)
library(emmeans)
library(purrr)
library(rlang)

# ---- Helper: build formula with interactions up to order k ----
build_formula <- function(y, facs, k) {
  facs <- unique(facs)
  # always include mains
  mains <- paste(facs, collapse = " + ")
  # interactions up to order k
  ints <- character(0)
  k <- min(k, length(facs))
  if (k >= 2) {
    for (ord in 2:k) {
      cmb <- combn(facs, ord, simplify = FALSE)
      ints <- c(ints, vapply(cmb, function(x) paste(x, collapse=":"), "", USE.NAMES = FALSE))
    }
  }
  rhs <- paste(c(mains, ints), collapse = " + ")
  as.formula(paste(y, "~", rhs))
}

# ---- Helper: drop factors with only 1 level (post-filter) ----
drop_constant_factors <- function(d, facs) {
  keep <- facs[vapply(facs, function(f) nlevels(factor(d[[f]])) > 1, logical(1))]
  unique(keep)
}

# ---- Main pipeline ----
run_anova_pipeline <- function(data,
                               response,              # e.g. "u_poc_h_1"
                               id_filter = NULL,      # e.g. "DMSPp" (optional)
                               time_filter = 24,      # or "t24" or NULL
                               factors = c("species","sal_adapt","sal_exp","time"),
                               max_interaction = 3,   # try up to 3-way, back off if needed
                               robust = TRUE,         # use White HC3 everywhere
                               posthoc_list = list(   # emmeans simple-effects specs
                                 list(factor = "sal_exp", by = c("species","sal_adapt"))
                               ),
                               verbose = TRUE,
                               transform = NULL       # function, e.g. log1p, or NULL
) {
  stopifnot(is.character(response), length(response) == 1)
  d <- data
  
  # -- Optional ID filter (e.g., keep DMSPp rows only)
  if (!is.null(id_filter) && "id" %in% names(d)) {
    d <- d %>% filter(id %in% id_filter)
  }
  
  # -- Optional time filter; accept 24 / "t24"
  if (!is.null(time_filter) && "time" %in% names(d)) {
    if (is.numeric(d$time) || is.integer(d$time)) {
      tf <- suppressWarnings(as.numeric(time_filter))
      d <- d %>% filter(time == tf)
    } else {
      d <- d %>% filter(time %in% time_filter)
    }
  }
  
  # -- Keep only rows with non-missing response and needed factors
  need_cols <- c(response, intersect(factors, names(d)))
  d <- d %>% select(any_of(unique(c(need_cols, "id"))))
  d <- d %>% filter(!is.na(.data[[response]]))
  
  # -- Coerce factors
  present_factors <- intersect(factors, names(d))
  for (f in present_factors) d[[f]] <- factor(d[[f]])
  
  # -- Drop factors that are constant after filtering
  present_factors <- drop_constant_factors(d, present_factors)
  
  if (length(present_factors) == 0) stop("No varying factors left after filtering.")
  
  # -- Optional transform
  y <- d[[response]]
  if (!is.null(transform)) {
    y <- transform(y)
  }
  d[[response]] <- y
  
  # -- Sum contrasts for Type III
  for (f in present_factors) {
    contrasts(d[[f]]) <- contr.sum(nlevels(d[[f]]))
  }
  
  # -- Try fitting from highest interaction order down until model is estimable
  k <- min(max_interaction, length(present_factors))
  fit <- NULL; formula_used <- NULL
  while (k >= 1 && is.null(fit)) {
    form <- build_formula(response, present_factors, k)
    fit_try <- try(lm(form, data = d), silent = TRUE)
    if (!inherits(fit_try, "try-error") && df.residual(fit_try) > 0) {
      fit <- fit_try
      formula_used <- form
    } else {
      k <- k - 1
    }
  }
  if (is.null(fit)) stop("Could not fit a model with positive residual df. Try fewer factors or more data.")
  
  if (verbose) {
    message("Fitted formula: ", deparse(formula_used))
    message("Residual df: ", df.residual(fit))
  }
  
  # -- Assumption checks (robustness will guard inference anyway)
  shapiro_p <- tryCatch({
    r <- residuals(fit); if (length(r) >= 3) shapiro.test(r)$p.value else NA_real_
  }, error = function(e) NA_real_)
  
  # Levene across all available groups (may fail if too many singletons)
  lev_all <- tryCatch({
    grp <- interaction(d[, present_factors], drop = TRUE, lex.order = TRUE)
    car::leveneTest(d[[response]], grp)
  }, error = function(e) NULL)
  
  # NCV (Breusch–Pagan)
  ncv_p <- tryCatch(car::ncvTest(fit)$p, error = function(e) NA_real_)
  
  # -- Type-III ANOVA (HC3 optional)
  if (robust) {
    aov_tab <- car::Anova(fit, type = 3, white.adjust = "hc3")
    vc <- car::hccm(fit, type = "hc3")
  } else {
    aov_tab <- car::Anova(fit, type = 3)
    vc <- vcov(fit)
  }
  
  # -- emmeans post-hocs (robust vcov carried into emmeans)
  posthoc <- list()
  if (length(posthoc_list)) {
    for (i in seq_along(posthoc_list)) {
      spec <- posthoc_list[[i]]
      fac <- spec$factor
      bys <- intersect(spec$by %||% character(0), present_factors)
      if (!(fac %in% present_factors)) next
      
      rhs <- if (length(bys)) paste0(fac, " | ", paste(bys, collapse = " * ")) else fac
      emm <- emmeans::emmeans(fit, specs = as.formula(paste("~", rhs)), vcov. = vc)
      # pairwise within the grid (Tukey adj)
      cmp <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
      posthoc[[rhs]] <- list(emm = emm, pairs = summary(cmp, infer = c(TRUE, TRUE)))
    }
  }
  
  list(
    data = d,
    response = response,
    factors_used = present_factors,
    formula = formula_used,
    fit = fit,
    anova_type3 = aov_tab,
    diagnostics = list(
      shapiro_p = shapiro_p,
      levene_all = lev_all,
      ncv_p = ncv_p
    ),
    posthoc = posthoc
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x


# Put this once near the top of the script
options(contrasts = c("contr.sum","contr.poly"))

run_anova2 <- function(dat, y, label) {
  stopifnot(is.character(y), length(y) == 1)
  # 1) Keep rows with non-missing outcome & required factors
  d <- dat %>%
    dplyr::filter(!is.na(.data[[y]])) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp)
  
  # 2) Check replicates per cell
  cell_vars <- c("species","sal_adapt","sal_exp")
  cell_n <- d %>%
    dplyr::count(dplyr::across(dplyr::all_of(cell_vars)), name = "n")
  
  # Decide whether to use means: singletons OR too few total rows
  need_means <- any(cell_n$n < 2) || nrow(d) < 8
  
  if (need_means) {
    bad <- cell_n %>% dplyr::filter(n < 2)
    if (nrow(bad) > 0) {
      msg_cells <- bad %>%
        dplyr::mutate(cell = paste(species, sal_adapt, sal_exp, sep="/")) %>%
        dplyr::pull(cell) %>% paste(collapse = ", ")
      message("Using cell means for ", label,
              " due to singleton cells (n<2): ", msg_cells)
    } else {
      message("Using cell means for ", label, " due to too few rows (n<8).")
    }
    
    d <- d %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cell_vars))) %>%
      dplyr::summarise(!!y := mean(.data[[y]], na.rm = TRUE), .groups = "drop")
  }
  
  if (nrow(d) < 8) {
    message("Too few rows for ", label, " even after checking/averaging.")
    return(NULL)
  }
  
  # 3) Fit the 3-way model
  form <- stats::as.formula(paste(y, "~ species * sal_adapt * sal_exp"))
  fit  <- stats::lm(form, data = d)
  
  # 4) Assumption checks
  resids <- residuals(fit)
  if (length(resids) < 3 || stats::sd(resids, na.rm = TRUE) < .Machine$double.eps) {
    shapiro_p <- NA_real_
    cat("Shapiro skipped: residuals have zero/near-zero variance or too few values.\n")
  } else if (length(resids) > 5000) {
    shapiro_p <- NA_real_
  } else {
    shapiro_p <- stats::shapiro.test(resids)$p.value
  }
  
  # Brown–Forsythe (median-centered) Levene tests
  lev3_p <- tryCatch(
    car::leveneTest(d[[y]] ~ interaction(d$species, d$sal_adapt, d$sal_exp), center = median)[[1,"Pr(>F)"]],
    error = function(e) NA_real_
  )
  lev1_p <- tryCatch(
    car::leveneTest(d[[y]] ~ d$sal_exp, center = median)[[1,"Pr(>F)"]],
    error = function(e) NA_real_
  )
  ncv_p  <- tryCatch(car::ncvTest(fit)$`p`, error = function(e) NA_real_)
  
  cat("\n============================\n", label, " — Assumption checks\n", sep="")
  cat("Shapiro (residuals) p = ", round(shapiro_p, 4), "\n", sep = "")
  cat("Levene (3-way groups, BF) p = ", round(lev3_p, 4), "\n", sep = "")
  cat("Levene (by sal_exp, BF)  p = ", round(lev1_p, 4), "\n", sep = "")
  cat("NCV test p = ", round(ncv_p, 4), "\n", sep = "")
  
  normal_ok <- is.na(shapiro_p) || shapiro_p > 0.05
  var_ok    <- (is.na(lev3_p) || lev3_p > 0.05) && (is.na(lev1_p) || lev1_p > 0.05)
  ncv_ok    <- is.na(ncv_p)    || ncv_p > 0.05
  use_hc3   <- !(normal_ok && var_ok && ncv_ok)
  
  cat(if (use_hc3) "\nUsing Type-III ANOVA with White HC3 robust SEs\n"
      else          "\nUsing standard Type-III ANOVA\n")
  
  aov_tab <- if (use_hc3) {
    tryCatch(car::Anova(fit, type = 3, white.adjust = "hc3"),
             error = function(e) {
               cat("\nHC3 failed (often leverage=1). Falling back to standard Type-III.\n")
               car::Anova(fit, type = 3)
             })
  } else {
    car::Anova(fit, type = 3)
  }
  print(aov_tab)
  
  # Warn about non-estimable terms (empty cells)
  ali <- tryCatch(stats::alias(fit), error = function(e) NULL)
  if (!is.null(ali) && (length(ali$Complete) + length(ali$Partial) > 0)) {
    message("Warning: Aliased (non-estimable) terms detected; interpret with care.")
  }
  
  # 5) Post-hoc Tukey: E25 vs E45 within species × sal_adapt, with the SAME robust vcov
  # ---- Post-hoc: E25 vs E45 within species × sal_adapt (robust vcov) ----
  cat("\nPost-hoc: sal_exp within species × sal_adapt\n")
  
  # robust covariance for emmeans (match HC3 used in omnibus)
  Vhc3 <- tryCatch(sandwich::vcovHC(fit, type = "HC3"),
                   error = function(e) tryCatch(car::hccm(fit, type = "hc3"),
                                                error = function(e) NULL))
  
  if (!is.null(Vhc3) && nlevels(droplevels(d$sal_exp)) > 1) {
    emm <- emmeans::emmeans(
      fit,
      ~ sal_exp | species * sal_adapt,
      vcov. = Vhc3
    )
    
    # EITHER of these is fine. Use ONE:
    
    # A) S3 method (no namespace)
    # cmp <- pairs(emm, adjust = "tukey")
    
    # B) Explicit exported function:
    cmp <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
    
    print(summary(cmp, infer = c(TRUE, TRUE), df = stats::df.residual(fit)))
  } else {
    message("Post-hoc not possible (no variation in sal_exp or robust vcov unavailable).")
  }
  
  
  invisible(list(fit = fit, anova = aov_tab))
}




# ---- Packages (add these) ----
library(dplyr)
library(lme4)        # mixed models
library(lmerTest)    # Satterthwaite/KR df for fixed-effect ANOVA tables
library(emmeans)

# ---- Set Type-III-safe contrasts once ----
options(contrasts = c("contr.sum","contr.poly"))

# ---- Mixed split-plot runner (replaces run_anova2) ----
run_anova2_mixed <- function(dat, y, label) {
  stopifnot(is.character(y), length(y) == 1)
  
  # Keep usable rows and ensure factors exist
  d <- dat %>%
    dplyr::filter(!is.na(.data[[y]])) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp) %>%
    mutate(
      species   = factor(species),
      sal_adapt = factor(sal_adapt),   # already A25/A45 upstream is fine too
      sal_exp   = factor(sal_exp),     # already E25/E45 upstream is fine too
      repl      = factor(repl),
      # bottle = biological replicate within species×A (whole plot)
      bottle    = interaction(species, sal_adapt, repl, drop = TRUE)
    )
  
  if (nrow(d) < 8L) {
    message("Too few rows for ", label, ".")
    return(invisible(NULL))
  }
  
  # Helpful sanity: do we have both E levels inside bottles?
  by_bottle <- d %>%
    group_by(bottle) %>%
    summarise(nE = n_distinct(droplevels(sal_exp)), .groups = "drop")
  nboth <- sum(by_bottle$nE >= 2)
  if (nboth == 0) {
    message("Skipping ", label, ": no bottle has both E levels; cannot estimate within-bottle E.")
    return(invisible(NULL))
  } else if (any(by_bottle$nE < 2)) {
    message("Note: some bottles miss one E level; E effect uses bottles with both levels where available.")
  }
  
  # ---- Mixed model: species * A * E + random intercept per bottle ----
  form <- stats::as.formula(paste(y, "~ species * sal_adapt * sal_exp + (1|species:sal_adapt:repl)"))
  fit  <- lmer(form, data = d, REML = TRUE)  # lmerTest provides F-tests with Satterthwaite df
  
  cat("\n============================\n", label, " — Mixed split-plot (REML)\n", sep = "")
  # Type-III ANOVA for fixed effects (valid with sum-to-zero contrasts)
  print(anova(fit, type = 3))   # use ddf="Kenward-Roger" if you prefer KR and have pbkrtest installed
  
  # ---- Post-hoc: E (E25 vs E45) within species × A (Tukey, mixed-model aware) ----
  cat("\nPost-hoc: E25 vs E45 within species × sal_adapt\n")
  emm <- emmeans(fit, ~ sal_exp | species * sal_adapt)
  cmp <- contrast(emm, method = "pairwise", adjust = "tukey")
  print(summary(cmp, infer = c(TRUE, TRUE)))  # lmerTest supplies df to emmeans
  
  invisible(list(fit = fit))
}



# ---- Packages ----
library(dplyr)
library(car)
library(emmeans)
library(lme4)
library(lmerTest)   # Satterthwaite/KR dfs for lmer
library(nlme)       # lme / gls with heteroscedasticity & correlation
library(e1071)      # skewness

options(contrasts = c("contr.sum","contr.poly"))

run_splitplot_adaptive <- function(dat, y, label, id_filter) {
  # Prep
  d0 <- dat %>%
    filter(time == "t24", id == id_filter, !is.na(.data[[y]])) %>%
    tidyr::drop_na(species, sal_adapt, sal_exp, repl) %>%
    mutate(
      species   = factor(species),
      sal_adapt = factor(sal_adapt),
      sal_exp   = factor(sal_exp),
      repl      = factor(repl),
      bottle    = interaction(species, sal_adapt, repl, drop = TRUE)
    )
  if (nrow(d0) < 8L) { message("Too few rows for ", label); return(invisible(NULL)) }
  
  # Keep bottles that inform within-bottle E effect (both E levels)
  d <- d0 %>% group_by(bottle) %>% filter(n_distinct(sal_exp) == 2) %>% ungroup()
  if (nrow(d) < 8L) d <- d0  # still allow model; post-hoc E may be limited
  
  notes <- c()
  transform <- NULL
  
  # --- quick data heuristics
  yvec <- d[[y]]
  pos_only <- all(yvec >= 0, na.rm = TRUE)
  sk <- tryCatch(e1071::skewness(yvec, na.rm = TRUE), error = function(e) NA_real_)
  
  # --- Helper: residual diag for a fitted model (works for lmer/lme/gls)
  diag_resid <- function(fit, data_for_groups) {
    r <- residuals(fit, type = "pearson")
    sh <- if (length(r) >= 3 && sd(r, na.rm = TRUE) > .Machine$double.eps)
      tryCatch(shapiro.test(as.numeric(r))$p.value, error = function(e) NA_real_)
    else NA_real_
    # Brown–Forsythe Levene on residuals by 3-way groups
    gr <- with(data_for_groups, interaction(species, sal_adapt, sal_exp, drop = TRUE))
    lev <- tryCatch(car::leveneTest(r ~ gr, center = median)[[1,"Pr(>F)"]], error = function(e) NA_real_)
    list(shapiro_p = sh, levene_p = lev)
  }
  
  # --- STEP 1: baseline LMM
  form_lmm <- as.formula(paste(y, "~ species * sal_adapt * sal_exp + (1|species:sal_adapt:repl)"))
  fit <- tryCatch(lmer(form_lmm, data = d, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                  error = function(e) NULL)
  singular <- tryCatch(isSingular(fit), error = function(e) TRUE)
  diag1 <- if (!is.null(fit)) diag_resid(fit, d) else list(shapiro_p = NA, levene_p = NA)
  
  # --- decide on transform (heavy skew + positive)
  if (!is.null(fit) && is.finite(sk) && sk > 1 & pos_only) {
    d[[paste0(y, "_log1p")]] <- log1p(d[[y]])
    y_tr <- paste0(y, "_log1p")
    fit2 <- tryCatch(lmer(as.formula(paste(y_tr, "~ species * sal_adapt * sal_exp + (1|species:sal_adapt:repl)")),
                          data = d, REML = TRUE,
                          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                     error = function(e) NULL)
    if (!is.null(fit2)) {
      fit <- fit2; transform <- "log1p"
      notes <- c(notes, "Applied log1p transform due to strong positive skew.")
      diag1 <- diag_resid(fit, d)
    }
  }
  
  # --- heteroscedasticity -> try LME with varIdent
  het_flag <- !is.na(diag1$levene_p) && diag1$levene_p < 0.05
  if (!is.null(fit) && (het_flag || singular)) {
    y_use <- if (is.null(transform)) y else paste0(y, "_log1p")
    form_fixed <- as.formula(paste(y_use, "~ species * sal_adapt * sal_exp"))
    # random = bottle nested in species:A
    lme_try <- function(weights_struct) {
      tryCatch(
        lme(fixed = form_fixed,
            random = ~1 | species/sal_adapt/repl,
            weights = weights_struct,
            data = d, method = "REML", control = lmeControl(msMaxIter = 200, msVerbose = FALSE)),
        error = function(e) NULL)
    }
    cand <- list(
      none      = lme_try(NULL),
      byE       = lme_try(varIdent(form = ~1 | sal_exp)),
      bySpecies = lme_try(varIdent(form = ~1 | species)),
      byA       = lme_try(varIdent(form = ~1 | sal_adapt)),
      byES      = lme_try(varComb(varIdent(form = ~1 | sal_exp), varIdent(form = ~1 | species)))
    )
    # choose best by AIC among successful fits
    ok <- cand[!vapply(cand, is.null, logical(1))]
    if (length(ok)) {
      aics <- sapply(ok, AIC)
      best_name <- names(which.min(aics))
      fit_lme <- ok[[best_name]]
      fit <- fit_lme
      notes <- c(notes, paste0("Switched to nlme::lme with variance structure: ", best_name, "."))
    }
  }
  
  # --- if still problematic (e.g., random effect collapses), use GLS with within-bottle correlation
  if (inherits(fit, "merMod")) {
    # we're fine; keep LMM
  } else if (inherits(fit, "lme")) {
    # fine; keep LME
  } else {
    # last resort GLS
    y_use <- if (is.null(transform)) y else paste0(y, "_log1p")
    fit_gls <- tryCatch(
      gls(as.formula(paste(y_use, "~ species * sal_adapt * sal_exp")),
          data = d,
          correlation = corCompSymm(form = ~1 | bottle),   # within-bottle correlation
          weights = varIdent(form = ~1 | sal_exp),         # allow diff σ by E
          method = "REML"),
      error = function(e) NULL)
    if (!is.null(fit_gls)) {
      fit <- fit_gls
      notes <- c(notes, "Fell back to gls with corCompSymm (within-bottle) + varIdent by E.")
    } else {
      message("Modeling failed for ", label)
      return(invisible(NULL))
    }
  }
  
  # --- Output: ANOVA for fixed effects
  cat("\n============================\n", label, " — Adaptive split-plot\n", sep = "")
  if (inherits(fit, "merMod")) {
    if (requireNamespace("pbkrtest", quietly = TRUE))
      print(anova(fit, type = 3, ddf = "Kenward-Roger"))
    else
      print(anova(fit, type = 3))  # Satterthwaite
  } else if (inherits(fit, "lme")) {
    print(anova(fit))  # Wald tests under the chosen variance structure
  } else if (inherits(fit, "gls")) {
    print(anova(fit))  # Wald tests for GLS
  }
  
  # --- Post-hoc: E within species × A
  cat("\nPost-hoc: E25 vs E45 within species × sal_adapt\n")
  emm <- emmeans(fit, ~ sal_exp | species * sal_adapt)
  print(summary(contrast(emm, "pairwise", adjust = "tukey"), infer = c(TRUE, TRUE)))
  
  # --- Decision notes
  if (!is.null(transform)) cat("\nDECISION:", transform, "transform used.\n")
  if (length(notes)) cat("DECISION:", paste(notes, collapse = "  "), "\n")
  
  invisible(list(fit = fit, notes = notes))
}
