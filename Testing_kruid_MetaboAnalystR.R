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



fls <- dir(path = "mzxml", full.names = TRUE)
pd <- data.frame(file = basename(fls),
                 sample = c("Peper_5", "Peper_6"),
                 group = "Peper")

peperdata <- readMSData(fls, pdata = new("NAnnotatedDataFrame", pd),
                   mode = "onDisk") 


MetaboAnalystR::ImportRawMSData(foldername = "mzxml", mode = "onDisk")

peperdata2 <- ImportRawMSData("C://Users/Micah/Documents/Data_IBD/", mode = 'onDisk', plotSettings = 'all')



oke <- PerformPeakProfiling(rawData = peperdata, Params = param_initial, plotSettings = SetPlotParam(Plot = TRUE))

kruiden_trimmed <- PerformDataTrimming("mzxml/KRUID_131/", rt.idx = 1)
param_initial <- SetPeakParam(platform = "general")
param_optimized <- PerformParamsOptimization(raw_data = kruiden_trimmed, param = param_initial)

parameter_data <- PerformDataTrimming('~/QC_IBD', rt.idx = '0.2')
param_initial <- SetPeakParam(platform = "UPLC-Q/E")
param_optimized <- PerformParamsOptimization(raw_data = parameter_data, param = param_initial)

psettings <- SetPlotParam(Plot = T,labels = T, format = 'png', dpi = 72)
ruwe__kruid_data <- ImportRawMSData(foldername = "mzxml", mode = "onDisk", plotSettings = psettings)
mSet <- PerformPeakProfiling(rawData = kruiden, Params = param_optimized)
annParams <- SetAnnotationParam(polarity = "positive")
annotPeaks <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams)
maPeaks <- FormatPeakList(annotPeaks = annotPeaks, annParams = annParams, filtIso = F, filtAdducts = F, missPercent = 1)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "metaboanalyst_input.csv", "colu", "disc")







