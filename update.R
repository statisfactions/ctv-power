library(dplyr)
library(stringr)
source("helper.R")

# Check available packages -----
available = rownames(available.packages())
pkgs = extract_cran("ctv-power.md")
not_there <- !(pkgs %in% available)

missing_pkgs = sort(unique(pkgs[not_there]))
missing_pkgs

# Correct any issues with cases -----
# case_replacements = available[tolower(available) %in% tolower(case_replacements)]
# names(case_replacements) = missing_pkgs[tolower(missing_pkgs) %in% tolower(available)]
# 
# correct_packages("ctv-power.md", case_replacements)

# render XML from md -------
out = md2ctv("ctv-power.md", links_file = "ctv-power-links.md")

# check compliance & render HTML ----------
ctv::ctv2html("PowerAnalysis.ctv")
