local({library(BiocManager); r <- repositories(); r["CRAN"] <- "http://cran.rstudio.com/"; options(repos = r)})
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                 paste(getRversion(), R.version$platform,
                       R.version$arch, R.version$os)))
# .libPaths("/usr/local/lib/R/library")