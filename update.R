source("helper.R")

# get contacts for feedback
pkgs = extract_pkgs("CRAN_task_view_fda.md")
not_there <- !(pkgs %in% rownames(available.packages()))
if (any(not_there)) {
  cat("###########################")
  cat(pkgs[not_there])
  cat("###########################")
}
to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_install)) install.packages(to_install)

cat(maintainers = unique(sapply(pkgs, get_maintainer_mail)))

# render XML from md
links <- c("<a href=\"http://www.psych.mcgill.ca/misc/fda/\">Website of the canonical FDA book</a>",
  "<a href=\"http://www.stat.ucdavis.edu/PACE/\">PACE: collection of MATLAB scripts from UC Davis</a>")
md2ctv("CRAN_task_view_fda.md", links)

# check compliance & render HTML
ctv::ctv2html("FDA.xml")
