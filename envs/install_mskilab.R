options(repos = c(CRAN = "https://cloud.r-project.org"))
for (pkg in c("mskilab/gUtils", "mskilab/bamUtils", "mskilab-org/gTrack", "mskilab/gGnome")) {
  remotes::install_github(pkg, upgrade = "never", dependencies = NA)
}
stopifnot(all(c("gUtils", "gTrack", "gGnome") %in% rownames(installed.packages())))
