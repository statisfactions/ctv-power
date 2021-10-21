library(dplyr)
library(stringr)
source("helper.R")

# get contacts for feedback

# Check available packages -----
available = rownames(available.packages())
pkgs = extract_pkgs("ctv-power.md")
not_there <- !(pkgs %in% available)

pkgs_tbl = data.frame(pkgs) %>% 
  count(pkgs) %>% 
  filter(n > 1) %>% 
  arrange(n) 

core.pkgs = 

missing_pkgs = sort(unique(pkgs[not_there]))
missing_pkgs

# Correct any issues with cases -----
# case_replacements = available[tolower(available) %in% tolower(case_replacements)]
# names(case_replacements) = missing_pkgs[tolower(missing_pkgs) %in% tolower(available)]
# 
# correct_packages("ctv-power.md", case_replacements)
# 
# 

if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}; BiocManager::install("ssizeRNA")

## CRAN installations

# to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
# 
# 
# if (length(to_install)) install.packages(to_install)
# 
# cat(maintainers = unique(sapply(pkgs, get_maintainer_mail)))

# render XML from md
links <- c("<a href=\"http://www.psych.mcgill.ca/misc/fda/\">Website of the canonical FDA book</a>",
  "<a href=\"http://www.stat.ucdavis.edu/PACE/\">PACE: collection of MATLAB scripts from UC Davis</a>")
md2ctv("CRAN_task_view_fda.md", links)

# check compliance & render HTML
ctv::ctv2html("FDA.xml")
