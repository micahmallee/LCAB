# Compound linking
## Data loading, database, peakprofiling etc
### LCAB:
data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
param_optimized <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0.01,
                           snthresh = 10, bw = 2, ppm = 20, min_peakwidth = 5, max_peakwidth = 30, 
                           noise = 100, prefilter = 3, value_of_prefilter = 100, minFraction = 0.5, 
                           minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                           family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                           mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized))
annParams <- SetAnnotationParam(polarity = 'positive')
annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams)

mona_msp <- read.msp(file = 'MoNA-export-GC-MS_Spectra.msp')

# rm(annParams, param_optimized, data)

## Create pseudospectra
xsetxcms <- xsAnnotate(xs = smSet$xcmsSet)
xsetxcms <- groupFWHM(object = xsetxcms)
### Correct format peakinfo
convertedgroupedxset <- xsetxcms
convertedgroupedxset@groupInfo <- convertedgroupedxset@xcmsSet@peaks
xsetxcmsmsp <- to.msp(object = convertedgroupedxset, file = NULL, settings = metaSetting(TSQXLS.GC, 'DBconstruction'))


## Compare pseudospectra to MoNA pseudospectra
allsimilarities <-  lapply(1:length(xsetmsp), function(x) {
  for (i in 1:length(DB.treated2)) {
    SpectrumSimilarity(spec.top = xsetmsp[[x]][, 1:2], spec.bottom = DB.treated2[[i]]$pspectrum, print.graphic = F, print.alignment = F)
  }
})


{ ## a peak table
  minI <- 0.1 * max(xsetxcms[1][, intensity])
  allpks <- xsetxcms[xsetxcms[, intensity] >= minI,]
  pspectra <- split(seq_len(nrow(allpks)), allpks[,"rt"])
}