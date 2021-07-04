if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
install.packages("pacman")
library(pacman)
library(devtools)

pacman::p_load("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", 
               "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph",
               "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR",
               "fgsea","httr","qs")
install.packages(c('glasso', 'crmn','plotly', 'magrittr', 'shinyFiles', 'shinycssloaders',
                   'shiny', 'shinyjs', 'stringr', 'shinyBS', 'DT',
                   'OrgMassSpecR', 'shinydashboard', 'crosstalk', 'shinydashboardPlus'))

BiocManager::install('metaMS')


devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
devtools::install_github("berlinguyinca/spectra-hash", subdir="splashR")
