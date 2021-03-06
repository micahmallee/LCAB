library(MetaboAnalystR)
library(xcms)
library(magrittr)

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
kruiden2 <- PerformDataTrimming(datapath = c('mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML',
                                  'mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), rt.idx = 1)

# Standaard parameters vaststellen
param_initial <- SetPeakParam(platform = "general", snthresh = 10)

# Door middel van trimmed kruiden data de parameters optimaliseren
param_optimized <- PerformParamsOptimization(raw_data = kruiden, param = param_initial)
param_optimized <- PerformParamsOptimization(raw_data = raw_kruiden_grouped, param = param_initial)

# Raw kruiden data inlezen
## Grouped data
raw_kruiden_grouped <- ImportRawMSData(foldername = 'mzxml/', mode = "onDisk", plotSettings = SetPlotParam(Plot = F))

## Solo data

# Peak profiling uitvoeren 
#'The PerformPeakProfiling function is an updated peak processing pipeline from XCMS R functions that performs peak detection, alignment, and grouping in an automatical step. 
#'The function also generates two diagnostic plots including statistics on the total intensity of peaks in different samples, a retention time adjustment map, 
#'and a PCA plot showing the overall sample clustering prior to data cleaning and statistical analysis.
mSet <- PerformPeakProfiling(rawData = raw_kruiden_grouped, Params = param_optimized2)


# Annotatie parameters vaststellen
annParams <- SetAnnotationParam(polarity = "positive")

# Peaklist maken die ingelezen kan worden door de Metaboanalyst webapp
annotPeaks2 <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams)
maPeaks <- FormatPeakList(annotPeaks = annotPeaks2, annParams = annParams, filtIso = F, filtAdducts = F, missPercent = 0, includeRT = T)


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
raw_peper <- readMSData(files = 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', mode = 'onDisk')
param_optimized <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0, snthresh = 10, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, noise = 0, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, mzCenterFun = "wMean")
# load('optimized_params')
param_optimized <- PerformParamsOptimization(raw_data = raw_peper, param = param_initial, ncore = 6)
smSet <- PerformPeakPicking(raw_peper, param = updateRawSpectraParam(param_optimized))
# smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized))
smSet <- MetaboAnalystR:::PerformPeakGrouping(smSet, param = updateRawSpectraParam(param_optimized))
# smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized))
annotPeaks_solo <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams2)

# Two samples:
raw_pepers <- readMSData(files = c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'kruiden/Kruid 132 Zwarte peper 7 191119me_72.mzXML'), mode = 'onDisk')
param_optimized2 <- PerformParamsOptimization(raw_data = raw_pepers, param = param_initial, ncore = 8)
smSet <- PerformPeakPicking(raw_pepers, param = updateRawSpectraParam(param_optimized2))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized2))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized2))
annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams2)

param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                snthresh = 10, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                noise = 1000, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                mzCenterFun = "wMean")

# Create function which makes a plotable xchromatogram object per sample with their respective found peaks
create_xchr <- function(mSet) {
  chr <- chromatogram(mSet[['onDiskData']])
  xchr <- as(chr, 'XChromatograms')
  chrompks <- mSet[["msFeatureData"]][["chromPeaks"]]
  chrompkd <- mSet[["msFeatureData"]][["chromPeakData"]]
  samples <- factor(chrompks[, "sample"], levels = 1:length(fileNames(mSet$onDiskData)))
  chrompks <- split.data.frame(chrompks, smpls)
  chrompkd <- split.data.frame(chrompkd, smpls)
  for (i in range(length(xchr))) {
    xchr[[i]]@chromPeaks <- chrompks[[i]]
    xchr[[i]]@chromPeakData <- chrompkd[[i]]
  }
  return(xchr)
}
