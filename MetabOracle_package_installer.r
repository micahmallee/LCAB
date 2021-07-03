if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")


install.packages("pacman")
library(pacman)
library(devtools)

pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", 
                 "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", 
                 "multtest","RBGL","edgeR","fgsea","httr","qs", "plotly", "shinyFiles", "xcms", "shinycssloaders",
                 "magrittr", "Rcpp", "CAMERA", "shiny", "shinyjs", "stringr", "mzR", "metaMS", "shinyBS", "DT", 
                 "OrgMassSpecR", "shinydashboard", "crosstalk", "digest", "ggplot2", "shinydashboardPlus"))




devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
devtools::install_github("berlinguyinca/spectra-hash", subdir="splashR")