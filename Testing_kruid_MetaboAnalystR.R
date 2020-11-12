library(MetaboAnalystR)
library(xcms)

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

# Eerst data trimmen/inlezen om daarmee later de parameters te bepalen (samples in 1 folder)
kruiden <- PerformDataTrimming("mzxml/", rt.idx = 1)

# Standaard parameters vaststellen
param_initial <- SetPeakParam(platform = "general", snthresh = 10)

# Door middel van trimmed kruiden data de parameters optimaliseren
param_optimized <- PerformParamsOptimization(raw_data = kruiden, param = param_initial)

# Raw kruiden data inlezen
## Grouped data
raw_kruiden <- ImportRawMSData(foldername = "mzxml/", mode = "onDisk", plotSettings = SetPlotParam(Plot = F))

## Solo data

# Peak profiling uitvoeren 
#'The PerformPeakProfiling function is an updated peak processing pipeline from XCMS R functions that performs peak detection, alignment, and grouping in an automatical step. 
#'The function also generates two diagnostic plots including statistics on the total intensity of peaks in different samples, a retention time adjustment map, 
#'and a PCA plot showing the overall sample clustering prior to data cleaning and statistical analysis.
mSet <- PerformPeakProfiling(rawData = raw_kruiden, Params = param_optimized)


# Annotatie parameters vaststellen
annParams <- SetAnnotationParam(polarity = "positive")

# Peaklist maken die ingelezen kan worden door de Metaboanalyst webapp
annotPeaks <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams)

# matched filter:
param_initial2 <- SetPeakParam(platform = 'general', Peak_method = 'matchedFilter')
param_optimized2 <- PerformParamsOptimization(raw_data = kruiden, param = param_initial2, ncore = 1)
#mSet2_matchedFilter <- PerformPeakProfiling(rawData = raw_kruiden, Params = param_optimized2$best_parameters)
mSet2_matchedFilter <- PerformPeakProfiling(rawData = raw_kruiden, Params = param_initial2, plotSettings = SetPlotParam(Plot = TRUE))
annParams2 <- SetAnnotationParam(polarity = "positive")
annotPeaks2 <- PerformPeakAnnotation(mSet = mSet2_matchedFilter, annotaParam = annParams2)

# cocosnoot:
cocosnoot_trimmed <- PerformDataTrimming("cocosnoot/", rt.idx = 1)
param_initial <- SetPeakParam(platform = "general", snthresh = 10)
param_optimized_cocosnoot <- PerformParamsOptimization(raw_data = cocosnoot_trimmed, param = param_initial)
raw_cocosnoot <- readMSData(files = 'cocosnoot/Kruid 46 Klapper_119.mzXML', mode = 'onDisk')

# one sample:
