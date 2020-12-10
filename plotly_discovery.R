# Plotly discovery
install.packages('plotly')
library(plotly)
library(MetaboAnalystR)
library(magrittr)
library(xcms)

data <- readMSData('mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML', mode = 'onDisk')
param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                 snthresh = 10, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                 noise = 0, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                 minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                 family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                 mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized2))

