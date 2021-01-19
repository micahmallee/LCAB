# Plotly discovery
install.packages('plotly')
library(plotly)
library(MetaboAnalystR)
library(magrittr)
library(xcms)
library(DT)
library(crosstalk)



# Osso:
data <- readMSData(c('mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML'), mode = 'onDisk')

# LCAB:
data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
data2 <- readMSData('kruiden/Kruid 130 Zwarte peper 5 191119me_70.mzXML', mode = 'onDisk')
param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                 snthresh = 100, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                 noise = 1000, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                 minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                 family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                 mzCenterFun = "wMean")

param_test <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0.01,
                             snthresh = 100, bw = 2, ppm = 22, min_peakwidth = 5, max_peakwidth = 30, 
                             noise = 1000, prefilter = 3, value_of_prefilter = 1000, minFraction = 0.5, 
                             minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                             family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                             mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data2, param = updateRawSpectraParam(param_test))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_test))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_test))


smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized2))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized2))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized2))

xchr <- create_xchr(smSet)
xchr1 <- create_xchr(smSet)

create_xchr <- function(mSet) {
  chr <- chromatogram(mSet[['onDiskData']])
  xchr <- as(chr, 'XChromatograms')
  chrompks <- mSet[["msFeatureData"]][["chromPeaks"]]
  chrompkd <- mSet[["msFeatureData"]][["chromPeakData"]]
  rt <- mSet[["msFeatureData"]][["adjustedRT"]]
  samples <- factor(chrompks[, "sample"], levels = 1:length(fileNames(mSet$onDiskData)))
  chrompks <- split.data.frame(chrompks, samples)
  chrompkd <- split.data.frame(chrompkd, samples)
  if (length(xchr) > 1) {
    for (i in 1:length(xchr)) {
      xchr[[i]]@rtime <- rt[[i]]
      xchr[[i]]@chromPeaks <- chrompks[[i]]
      xchr[[i]]@chromPeakData <- chrompkd[[i]]
    }
  } else {
    for (i in 1:length(xchr)) {
      xchr[[i]]@chromPeaks <- chrompks[[i]]
      xchr[[i]]@chromPeakData <- chrompkd[[i]]
    }
  }
  return(xchr)
}

p <- plot_ly(x =  xchr[[1]]@rtime, y = xchr[[1]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities', source = 'peakplot')
p <- p %>% add_trace(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], 
                     type = 'scatter', mode = 'markers', name = paste0('Sample ', 1), text = xchr[[1]]@chromPeaks[,1], 
                     hoverinfo = 'text') %>% highlight('plotly_selected', dynamic = F)
p




# BS
fig <- plot_ly(x =  xchr[[1]]@rtime, y = xchr[[1]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities')
fig <- fig %>% add_trace(x =  xchr[[2]]@rtime, y = xchr[[2]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities2')
fig <- fig %>% add_trace(x =  xchr[[2]]@rtime, y = xchr[[2]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities2')
fig <- fig %>% add_trace(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], type = 'contour', mode = 'markers', name = 'Sample 1')
fig <- fig %>% add_trace(x =  xchr[[2]]@chromPeaks[,4], y = xchr[[2]]@chromPeaks[,9], type = 'scatter', mode = 'markers', name = 'Sample 2')
fig

fig <- plot_ly(x =  xchr[[1]]@rtime, y = xchr[[1]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities')
fig <- plot_ly(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], type = 'bar', name = 'Sample 1')
fig
rm(fig)

fig <- plot_ly(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], type = 'scatter', mode = 'markers+lines')

fig <- plot_ly(x =  smSet$onDiskData@featureData@data$retentionTime, y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt 1')
fig <- fig %>% add_trace(x =  smSet$xcmsSet@rt$raw, y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt raw2')
fig <- fig %>% add_trace(x =  smSet[["xcmsSet"]]@rt[["corrected"]][["1"]], y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt corrected1')

fig

# Oke wat moet er gebeuren?
#' chromatogram met tic of bpi, met daarop de gevonden pieken, per sample andere kleur.
#' Hover geeft meer informatie over de individuele pieken.
#' eventueel klikken op pieken het massa spectrum te zien?

# Starten met de base:
fig <- plot_ly(x =  smSet$onDiskData@featureData@data$retentionTime, y = unlist(spectrapply(smSet$onDiskData, FUN = function(z) sum(z@intensity))), mode = 'lines')
fig <- plot_ly(x =  smSet$onDiskData@featureData@data$retentionTime, y = smSet$onDiskData@featureData@data$totIonCurrent, mode = 'lines')

# oke <- filterRt(smSet$, rt = c(570, 600))



fig
fig <- plot_ly(x = smSet$msFeatureData$chromPeaks[,4], y = smSet$msFeatureData$chromPeaks[,9], type = 'scatter', mode = 'markers', name = 'foundpeaks', 
               text = smSet$msFeatureData$chromPeaks[,1], hoverinfo = 'text')
fig <- plot_ly(x = smSet$xcmsSet@peaks[,4], y = smSet$xcmsSet@peaks[,9], type = 'scatter', mode = 'markers', name = 'foundpeaks')

nr <- nrow(smSet$onDiskData)
nc <- ncol(smSet$onDiskData)

if (nr > 1) {
  par(mfrow = c(round(sqrt(nr)), ceiling(sqrt(nr))))
}
for (i in seq_len(nr)) {
  if (nc > 1) 
    .plotChromatogramList(smSet$onDiskData[i, , drop = TRUE])
  else plot(smSet$onDiskData[i, 1])
}


oke <- filterRt(smSet$onDiskData, rt = c(1824.569, 1824.571)) %>% filterMz(mz = c(122.06517, 122.06718))
plot(oke)

p <- ggplot(data.frame(x = rt, aes(x = x, ymax = y, ymin = 0)) + geom_linerange() + ylab(label = 'Intensity') + xlab('Retention time'))
plotteke <- ggplot(data.frame(x = mz, y = ints))
ggplotly(plotteke)
ggplotly(p)




