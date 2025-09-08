
library(readxl)
library(janitor)
library(readr)

#setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")

#RERUNS!! 
pol <- read_excel("data/raw/Data_Final_FMK.xlsx", sheet = 1, skip = 4)
chr <- read_excel("data/raw/Data_Final_FMK.xlsx", sheet = 2, skip = 4) 

#pol <- read_excel("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/clean/ye.xlsx", sheet = 1, skip = 4) 
#chr <- read_excel("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/clean/ye.xlsx", sheet = 2, skip = 4) 

################################################################
#cleaning up sheet

chr <- chr %>% janitor::clean_names(replace = c("\u00b5" = "u"))
pol <- pol %>% janitor::clean_names(replace = c("\u00b5" = "u"))



#polarella, sheet 1
#pol <- clean_names(pol)
#pol <- remove_empty(pol, which = c("rows", "cols"))
#write_csv(pol, "clean_pol.csv")

#chryso, sheet 2
#chr <- clean_names(chr)
#chr <- remove_empty(chr, which = c("rows", "cols"))
#write_csv(chr, "clean_chr.csv")


pol <- clean_names(pol)
pol <- remove_empty(pol, which = c("rows", "cols"))
write_csv(pol, file.path("data", "processed", "clean_pol.csv"))

# chryso, sheet 2
chr <- clean_names(chr)
chr <- remove_empty(chr, which = c("rows", "cols"))
write_csv(chr, file.path("data", "processed", "clean_chr.csv"))

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
write_csv(pol, file.path("data", "processed", "clean_pol_modified.csv"))


###############################################################

#View(pol)
#View(chr)
