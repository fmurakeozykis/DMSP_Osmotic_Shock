# ==========================================
# Uptake-family figures (raw t24 boxplots)
# D3-DMSP taken up  |  D3-DMS(P) lost/demethylated
# ==========================================
setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")
library(dplyr)
library(readr)
library(ggplot2)
library(emmeans)   # emmeans::cld for post-hoc letters

# ---------- LOAD + PREP (skip if df already exists) ----------
# chr <- read_csv("clean_chr.csv")
# pol <- read_csv("clean_pol_modified.csv")  # or clean_pol.csv
# df  <- bind_rows(chr, pol)

# Keep needed vars & standardise time in one step
# Comment this block if df already prepared
keep_cols <- c(
  "species","sal_adapt","sal_exp","time","repl",
  "d3_p_taken_up","d3_p_lost_demethylated"
)
df <- df %>%
  dplyr::select(dplyr::any_of(keep_cols)) %>%
  mutate(
    time = factor(
      ifelse(as.character(time) %in% c("0","t0"), "t0",
             ifelse(as.character(time) %in% c("24","t24"), "t24", NA_character_)),
      levels = c("t0","t24")
    ),
    sal_adapt = factor(sal_adapt),
    sal_exp   = factor(sal_exp),
    group_x   = interaction(species, sal_adapt, sep = " × ")
  ) %>%
  filter(!is.na(sal_adapt), !is.na(sal_exp))
# -------------------------------------------------------------

# ---------- Helper: t24 boxplot + (optional) Holm letters ----------
make_uptake_plot <- function(dat, response_var, y_lab) {
  dat_t24 <- dat %>%
    filter(time == "t24") %>%
    mutate(group_x = interaction(species, sal_adapt, sep = " × "))
  
  p <- ggplot(dat_t24, aes(x = group_x, y = .data[[response_var]], fill = sal_exp)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(aes(color = sal_exp), width = 0.15, alpha = 0.6,
                size = 1.6, show.legend = FALSE) +
    labs(x = "Species × acclimation salinity",
         y = y_lab,
         fill = "Exposure salinity") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),
          legend.position = "top")
  
  # ----- OPTIONAL: Holm-adjusted post-hoc letters (comment to disable) -----
  try({
    mod <- lm(as.formula(paste(response_var, "~ species * sal_adapt * sal_exp")),
              data = dat_t24)
    emm <- emmeans(mod, ~ species * sal_adapt * sal_exp)
    cld_df <- emmeans::cld(emm, adjust = "holm", Letters = letters) |>
      as.data.frame() |>
      mutate(group_x = interaction(species, sal_adapt, sep = " × "))
    
    # place letters just above local max per x-group
    y_pos <- dat_t24 |>
      group_by(group_x) |>
      summarise(y_top = max(.data[[response_var]], na.rm = TRUE) * 1.05,
                .groups = "drop")
    
    cld_df <- cld_df |>
      left_join(y_pos, by = "group_x")
    
    p <- p + geom_text(data = cld_df,
                       aes(x = group_x, y = y_top, label = .group),
                       inherit.aes = FALSE)
  }, silent = TRUE)
  # ------------------------------------------------------------------------
  
  p
}

# ---------- Make the two uptake plots ----------
p_taken   <- make_uptake_plot(df, "d3_p_taken_up",          "D3-DMSP taken up (raw units)")
p_lostdem <- make_uptake_plot(df, "d3_p_lost_demethylated", "D3-DMS(P) lost/demethylated (raw units)")

# View
print(p_taken)
print(p_lostdem)

# Save (optional)
# ggsave("fig_D3_taken_up_t24.png",          p_taken,   width = 7, height = 4.5, dpi = 300, bg = "white")
# ggsave("fig_D3_lost_demethylated_t24.png", p_lostdem, width = 7, height = 4.5, dpi = 300, bg = "white")





















# ==========================================
# Build clean dataframe for uptake plots
# ==========================================

library(dplyr)
library(readr)

# Load (uncomment if not loaded yet)
chr <- read_csv("clean_chr.csv")
pol <- read_csv("clean_pol_modified.csv")  # or clean_pol.csv

# Start fresh: don't use the name `df`
dat_all <- bind_rows(chr, pol)

keep_cols <- c(
  "species","sal_adapt","sal_exp","time","repl",
  "d3_p_taken_up","d3_p_lost_demethylated"
)

dat_uptake <- dat_all %>%
  dplyr::select(dplyr::any_of(keep_cols)) %>%
  dplyr::mutate(
    time = factor(
      dplyr::case_when(
        as.character(time) %in% c("0","t0")  ~ "t0",
        as.character(time) %in% c("24","t24")~ "t24",
        TRUE ~ NA_character_
      ),
      levels = c("t0","t24")
    ),
    sal_adapt = factor(sal_adapt),
    sal_exp   = factor(sal_exp),
    group_x   = interaction(species, sal_adapt, sep = " × ")
  ) %>%
  dplyr::filter(!is.na(sal_adapt), !is.na(sal_exp))

# sanity check
dplyr::glimpse(dat_uptake)


library(ggplot2)
library(emmeans)

make_uptake_plot <- function(dat, response_var, y_lab) {
  dat_t24 <- dat %>%
    dplyr::filter(time == "t24") %>%
    dplyr::mutate(group_x = interaction(species, sal_adapt, sep = " × "))
  
  p <- ggplot(dat_t24, aes(x = group_x, y = .data[[response_var]], fill = sal_exp)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(aes(color = sal_exp), width = 0.15, alpha = 0.6,
                size = 1.6, show.legend = FALSE) +
    labs(x = "Species × acclimation salinity", y = y_lab, fill = "Exposure salinity") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), legend.position = "top")
  
  # ----- OPTIONAL: post-hoc letters (comment to disable) -----
  try({
    mod <- lm(as.formula(paste(response_var, "~ species * sal_adapt * sal_exp")),
              data = dat_t24)
    emm <- emmeans(mod, ~ species * sal_adapt * sal_exp)
    cld_df <- emmeans::cld(emm, adjust = "holm", Letters = letters) |>
      as.data.frame() |>
      dplyr::mutate(group_x = interaction(species, sal_adapt, sep = " × "))
    
    y_pos <- dat_t24 |>
      dplyr::group_by(group_x) |>
      dplyr::summarise(y_top = max(.data[[response_var]], na.rm = TRUE) * 1.05,
                       .groups = "drop")
    
    cld_df <- dplyr::left_join(cld_df, y_pos, by = "group_x")
    
    p <- p + geom_text(data = cld_df,
                       aes(x = group_x, y = y_top, label = .group),
                       inherit.aes = FALSE)
  }, silent = TRUE)
  # -----------------------------------------------------------
  
  p
}

# Make plots
p_taken   <- make_uptake_plot(dat_uptake, "d3_p_taken_up",          "D3-DMSP taken up (% of total D3P added)")
p_lostdem <- make_uptake_plot(dat_uptake, "d3_p_lost_demethylated", "D3-DMS(P) lost/demethylated (% of total D3P added)")

print(p_taken)
print(p_lostdem)
















# ==========================================
# Build clean dataframe for uptake plots (DMS(P)t only)
# ==========================================

library(dplyr)
library(readr)

# Load (uncomment if not loaded yet)
# chr <- read_csv("clean_chr.csv")
# pol <- read_csv("clean_pol_modified.csv")  # or clean_pol.csv

# Start fresh: don't use the name `df` to avoid clashes
dat_all <- bind_rows(chr, pol)

keep_cols <- c(
  "species","sal_adapt","sal_exp","time","repl","id",
  "d3_p_taken_up","d3_p_lost_demethylated"
)

dat_uptake <- dat_all %>%
  dplyr::filter(id == "DMS(P)t") %>%
  dplyr::select(dplyr::any_of(keep_cols)) %>%
  dplyr::mutate(
    time = factor(
      dplyr::case_when(
        as.character(time) %in% c("0","t0")  ~ "t0",
        as.character(time) %in% c("24","t24")~ "t24",
        TRUE ~ NA_character_
      ),
      levels = c("t0","t24")
    ),
    sal_adapt = factor(sal_adapt),
    sal_exp   = factor(sal_exp),
    group_x   = interaction(species, sal_adapt, sep = " × ")
  ) %>%
  dplyr::filter(!is.na(sal_adapt), !is.na(sal_exp))

# sanity check
dplyr::glimpse(dat_uptake)

library(ggplot2)
library(emmeans)

make_uptake_plot <- function(dat, response_var, y_lab) {
  dat_t24 <- dat %>%
    dplyr::filter(time == "t24") %>%
    dplyr::mutate(group_x = interaction(species, sal_adapt, sep = " × "))
  
  p <- ggplot(dat_t24, aes(x = group_x, y = .data[[response_var]], fill = sal_exp)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(aes(color = sal_exp), width = 0.15, alpha = 0.6,
                size = 1.6, show.legend = FALSE) +
    labs(x = "Species × acclimation salinity", y = y_lab, fill = "Exposure salinity (ppt)") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), legend.position = "top")
  
  # ----- OPTIONAL: post-hoc letters (comment to disable) -----
  try({
    mod <- lm(as.formula(paste(response_var, "~ species * sal_adapt * sal_exp")),
              data = dat_t24)
    emm <- emmeans(mod, ~ species * sal_adapt * sal_exp)
    cld_df <- emmeans::cld(emm, adjust = "holm", Letters = letters) |>
      as.data.frame() |>
      dplyr::mutate(group_x = interaction(species, sal_adapt, sep = " × "))
    
    y_pos <- dat_t24 |>
      dplyr::group_by(group_x) |>
      dplyr::summarise(y_top = max(.data[[response_var]], na.rm = TRUE) * 1.05,
                       .groups = "drop")
    
    cld_df <- dplyr::left_join(cld_df, y_pos, by = "group_x")
    
    p <- p + geom_text(data = cld_df,
                       aes(x = group_x, y = y_top, label = .group),
                       inherit.aes = FALSE)
  }, silent = TRUE)
  # -----------------------------------------------------------
  
  p
}

# Make plots
p_taken   <- make_uptake_plot(dat_uptake, "d3_p_taken_up",          "D3-DMSP taken up (% of D3P added)")
p_lostdem <- make_uptake_plot(dat_uptake, "d3_p_lost_demethylated", "D3-DMS(P) lost/demethylated (% of D3P added)")

print(p_taken)
print(p_lostdem)


