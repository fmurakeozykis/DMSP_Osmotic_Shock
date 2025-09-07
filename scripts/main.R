
library(readxl)
library(janitor)
library(readr)

setwd("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/R/csvs_for_chat")

#RERUNS!!
pol <- read_excel("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/clean/ye.xlsx", sheet = 1, skip = 4) 
chr <- read_excel("C:/Users/fmura/Documents/groningen/Hon Project Arctic Algae/Data/clean/ye.xlsx", sheet = 2, skip = 4) 

################################################################
#cleaning up sheet

chr <- chr %>% janitor::clean_names(replace = c("\u00b5" = "u"))
pol <- pol %>% janitor::clean_names(replace = c("\u00b5" = "u"))



#polarella, sheet 1
pol <- clean_names(pol)
pol <- remove_empty(pol, which = c("rows", "cols"))
write_csv(pol, "clean_pol.csv")

#chryso, sheet 2
chr <- clean_names(chr)
chr <- remove_empty(chr, which = c("rows", "cols"))
write_csv(chr, "clean_chr.csv")

###############################################################

#View(pol)
#View(chr)
