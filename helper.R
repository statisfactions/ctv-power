## Input:
##  file  [character] file path to taskview .md file.
## Returns:
##  vector with package names that are listed in the task view.
extract_cran = function(file) {
  ll = readLines(file)
  pkg_strings1 = ll[grepl("\\[\\S*?\\]\\(https://CRAN", ll)]
  
  ## extract package names
  pkg_strings2 = unlist(stringr::str_extract_all(pkg_strings1, "\\[\\S*?\\]\\(https://CRAN"))
  pkgs = gsub(".*\\[(\\S*)\\].*", "\\1", pkg_strings2)
  
  ## extract links to CRAN
  #links = gsub(".*\\]\\((\\S*)\\).*", "\\1", pkg_strings)
  
  ## return package names
  pkgs
}

extract_github = function(file) {
  ll = readLines(file)
  pkg_strings1 = ll[grepl("\\[\\S*?\\]\\(https://github", ll)]
  
  ## extract package names
  pkg_strings2 = unlist(stringr::str_extract_all(pkg_strings1, "\\[\\S*?\\]\\(https://github"))
  pkgs = gsub(".*\\[(\\S*)\\].*", "\\1", pkg_strings2)
  
  ## extract links to CRAN
  #links = gsub(".*\\]\\((\\S*)\\).*", "\\1", pkg_strings)
  
  ## return package names
  pkgs
}

## Corrects mismatched packages, given a vector of replacements.
## Mainly used right now for corrected case issues
## Rewrites out to same file
correct_packages = function(file, replacements) {
  ll = readLines(file)
  
  ll_out = stringr::str_replace_all(ll, replacements)
  
  writeLines(ll_out, con = file)
  
}

## function that, given a package name, returns the email of the maintainer
## of that package
get_maintainer_mail = function(pkg_name) {
  desc = packageDescription(pkg_name)
  sub(".*<(.*)>.*", "\\1", desc$Maintainer)
}

# turn md into ctv-compliant XML
md2ctv <- function(file = "ctv-power.md", links_file = "ctv-power-links.md", save = TRUE) {
  
  md2htmlclean = function(file) {
    html = rmarkdown::render(file)
    ll = readLines(html)
    ## remove HTML header & footer
    headerfooter = grep("<!--", ll)
    ll = ll[(headerfooter[1] + 1) : (headerfooter[2] - 1)]
    ## remove comments
    ll = ll[-grep("-->", ll)]
    ## remove div
    div = grep("<[/]?div", ll)
    if(length(div))
      ll = ll[-div]
    ## replace headers...
    # ll = sub("<h4>((\\w+ ?)*).*</h4>$", replacement = "<p><strong>\\1</strong></p>", ll)
    ## ... and package hyperlinks
    ll = sub("<a href=\"https://CRAN.+>(\\w+)</a>", replacement = "<pkg>\\1</pkg>", ll)
  }
  
  info = md2htmlclean(file)
  links = md2htmlclean(links_file)
  pkgs = extract_pkgs(file)
  
  pkgcounts = as.matrix(table(pkgs))[,1]
  uniq_pkgs = sort(unique(pkgs))
  which_core = names(pkgcounts)[pkgcounts >= 5]
  
  pkgstrings = ifelse(uniq_pkgs %in% which_core,
         paste0('<pkg priority="core">', uniq_pkgs, "</pkg>"),
         paste0("<pkg>", uniq_pkgs, "</pkg>"))
  
  # add XML header & footer
  pkglist = c("<packagelist>", paste0("<pkg>", uniq_pkgs, "</pkg>"), "</packagelist>")
  ll = c("<CRANTaskView>", info, "</info>", pkglist,
         "<links>", links, "</links>",
         "</CRANTaskView>")
  
  writeLines(ll, "PowerAnalysis.ctv")
  
  invisible(ll)
}



