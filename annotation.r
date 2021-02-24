# Mummichog for marathon project

data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
# data2 <- readMSData('kruiden/Kruid 130 Zwarte peper 5 191119me_70.mzXML', mode = 'onDisk')
param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                 snthresh = 100, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                 noise = 1000, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                 minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                 family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                 mzCenterFun = "wMean")

# param_test <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0.01,
#                            snthresh = 10, bw = 2, ppm = 22, min_peakwidth = 5, max_peakwidth = 30, 
#                            noise = 1000, prefilter = 3, value_of_prefilter = 0.02, minFraction = 0.5, 
#                            minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
#                            family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
#                            mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized2))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized2))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized2))
smSet$xcmsSet@phenoData$sample_group <- gsub(' ', x = smSet$xcmsSet@phenoData$sample_name, replacement = '_')

annParams <- SetAnnotationParam(polarity = 'positive')
annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams)


maPeaks <- FormatPeakList(annotPeaks = annotPeaks, annParams = annParams, filtIso = F, filtAdducts = F, missPercent = 0.5, includeRT = T)


smSet <- InitDataObjects("pktable", 'stat', paired = FALSE)
smSet <- Read.TextData(mSetObj = smSet, filePath = 'metaboanalyst_input.csv', format = 'colu', lbl.type = 'disc')



