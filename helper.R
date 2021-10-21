## Input:
##  file  [character] file path to taskview .md file.
## Returns:
##  vector with package names that are listed in the task view.
extract_pkgs = function(file) {
  ll = readLines(file)
  pkg_strings = ll[grepl("\\[\\S*\\]\\(https://cran.*\\)", ll)]

  ## extract package names
  pkgs = gsub(".*\\[(\\S*)\\].*", "\\1", pkg_strings)

  ## extract links to CRAN
  #links = gsub(".*\\]\\((\\S*)\\).*", "\\1", pkg_strings)

  ## return package names
  pkgs
}

## function that, given a package name, returns the email of the maintainer
## of that package
get_maintainer_mail = function(pkg_name) {
  desc = packageDescription(pkg_name)
  sub(".*<(.*)>.*", "\\1", desc$Maintainer)
}

# turn md into ctv-compliant XML
md2ctv <- function(file, links, save = TRUE) {
  pkgs = extract_pkgs(file)
  html = rmarkdown::render(file)
  ll = readLines(html)
  ## remove HTML header & footer
  headerfooter = grep("<!--", ll)
  ll = ll[(headerfooter[1] + 1) : (headerfooter[2] - 1)]
  ## remove comments
  ll = ll[-grep("-->", ll)]
  ## remove div
  ll = ll[-grep("<[/]?div", ll)]
  ## replace headers...
  ll = sub("<h4>((\\w+ ?)*).*</h4>$", replacement = "<p><strong>\\1</strong></p>", ll)
  ## ... and package hyperlinks
  ll = sub("<a href=\"https://cran.+>(\\w+)</a>", replacement = "<pkg>\\1</pkg>", ll)

  # add XML header & footer
  pkglist = c("<packagelist>", paste0("<pkg>", sort(unique(pkgs)), "</pkg>"), "</packagelist>")
  ll = c("<CRANTaskView>", ll, "</info>", pkglist,
    "<links>", links, "</links>",
    "</CRANTaskView>")

  if (save) writeLines(ll, "FDA.xml")
  invisible(ll)
}



