# Plotly discovery
install.packages('plotly')
library(plotly)
library(MetaboAnalystR)
library(magrittr)
library(xcms)

data <- readMSData(c('mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML'), mode = 'onDisk')
param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                 snthresh = 100, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                 noise = 0, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                 minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                 family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                 mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized2))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized2))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized2))

xchr <- create_xchr(smSet)

chrmaken <- function (rtime = numeric(), intensity = numeric(), mz = c(NA_real_, 
                                                           NA_real_), filterMz = c(NA_real_, NA_real_), precursorMz = c(NA_real_, 
                                                                                                                        NA_real_), productMz = c(NA_real_, NA_real_), fromFile = integer(), 
          aggregationFun = character(), msLevel = 1L) 
{
  if (is.unsorted(rtime)) {
    idx <- order(rtime)
    rtime <- rtime[idx]
    intensity <- intensity[idx]
  }
  new("Chromatogram", rtime = rtime, intensity = intensity, 
      mz = range(mz), filterMz = range(filterMz), precursorMz = range(precursorMz), 
      productMz = range(productMz), fromFile = as.integer(fromFile), 
      aggregationFun = aggregationFun, msLevel = as.integer(msLevel))
}



a <- list(
  autotick = F
)

fig <- plot_ly()


create_xchr <- function(mSet) {
  chr <- chromatogram(mSet[['onDiskData']])
  xchr <- as(chr, 'XChromatograms')
  chrompks <- mSet[["msFeatureData"]][["chromPeaks"]]
  chrompkd <- mSet[["msFeatureData"]][["chromPeakData"]]
  samples <- factor(chrompks[, "sample"], levels = 1:length(fileNames(mSet$onDiskData)))
  chrompks <- split.data.frame(chrompks, samples)
  chrompkd <- split.data.frame(chrompkd, samples)
  for (i in 1:length(xchr)) {
    xchr[[i]]@chromPeaks <- chrompks[[i]]
    xchr[[i]]@chromPeakData <- chrompkd[[i]]
  }
  return(xchr)
}

a <- list(
  tick0 = 50
)

fig <- plot_ly(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], type = 'scatter', mode = 'markers+lines')
fig <- plot_ly(x =  xchr[[1]]@rtime, y = xchr[[1]]@intensity, type = 'scatter', mode = 'lines', name = 'base intensities')
fig <- plot_ly(x =  smSet$onDiskData@featureData@data$retentionTime, y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt 1')
fig <- fig %>% add_trace(x =  smSet$xcmsSet@rt$raw, y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt raw2')
fig <- fig %>% add_trace(x =  smSet[["xcmsSet"]]@rt[["corrected"]][["1"]], y = smSet$onDiskData@featureData@data$basePeakIntensity, type = 'scatter', mode = 'lines', name = 'rt corrected1')
fig <- fig %>% add_trace(x =  xchr[[1]]@chromPeaks[,4], y = xchr[[1]]@chromPeaks[,9], type = 'scatter', mode = 'markers+lines', name = 'Sample 1')
fig <- fig %>% add_trace(x =  xchr[[2]]@chromPeaks[,4], y = xchr[[2]]@chromPeaks[,9], type = 'scatter', mode = 'markers+lines', name = 'Sample 2')
fig

fig <- plot_ly(x = smSet$msFeatureData$chromPeaks[,4], y = smSet$msFeatureData$chromPeaks[,9], type = 'scatter', mode = 'markers', name = 'foundpeaks', text = smSet$msFeatureData$chromPeaks[,1], hoverinfo = 'text')
fig <- plot_ly(x = smSet$xcmsSet@peaks[,4], y = smSet$xcmsSet@peaks[,9], type = 'scatter', mode = 'markers', name = 'foundpeaks')

rm(fig)

oke <- filterRt(smSet$onDiskData, rt = c(618, 619))
plot(oke)


p <- ggplot(x = rt, y = ints)
ggplot(p)

# data frame nodig met rt, peaks, origin sample, mz, intensity
# basic info
rt <- rtime(smSet$onDiskData)
mz <- mz(smSet$onDiskData)
ints <- intensity(smSet$onDiskData)
# Peak info
pks <- smSet$msFeatureData$chromPeaks
pkd <- smSet$msFeatureData$chromPeakData
smpls <- factor(pks[, "sample"], levels = 1:length(fileNames(smSet$onDiskData)))
pks <- split.data.frame(pks, smpls)
pkd <- split.data.frame(pkd, smpls)

fig <- plot_ly(data = smSet$msFeatureData$chromPeaks, x = rt, y = ints, type = 'scatter', mode = 'markers+lines')


p <- ggplot(data.frame(x = rt, aes(x = x, ymax = y, ymin = 0)) + geom_linerange() + ylab(label = 'Intensity') + xlab('Retention time'))
plotteke <- ggplot(data.frame(x = mz, y = ints))
ggplotly(plotteke)
ggplotly(p)


pksall <- chromPeaks(xchr[[1]]) # grab all peaks
pksnr <- nrow(pksall) # count peaks
peakCol <- "#00000060"
peakBg <- "#00000020"
peakPch <- 1
peakCol <- rep(peakCol[1], pksnr) # give each peak a color
peakBg <- rep(peakBg[1], pksnr)
peakPch <- rep(peakPch[1], pksnr)


nr <- pksnr
x <- xchr[[1]]
for (i in seq_len(nr)) {
  x_sub <- x[i, , drop = FALSE]
  plot(as(x_sub, ifelse(is(x_sub, "XChromatograms"), 
                        "Chromatograms", "Chromatogram")), col = col, 
       lty = lty, type = type, xlab = xlab, ylab = ylab, 
       main = main, ...)
  idx <- which(pks_all[, "row"] == i)
  if (length(idx) && peakType != "none") {
    pks <- chromPeaks(x_sub)
    .add_chromatogram_peaks(x_sub, pks, col = peakCol[idx], 
                            bg = peakBg[idx], type = peakType, pch = peakPch[idx], 
                            ...)
  }
}




xchrplotten <- function (x, y, ...) 
{
  .local <- function (x, col = "#00000060", lty = 1, type = "l", 
                      xlab = "retention time", ylab = "intensity", main = NULL, 
                      peakType = c("polygon", "point", "rectangle", "none"), 
                      peakCol = "#00000060", peakBg = "#00000020", peakPch = 1, 
                      ...) 
  {
    peakType <- match.arg(peakType)
    nr <- nrow(x)
    if (nr > 1) 
      par(mfrow = c(round(sqrt(nr)), ceiling(sqrt(nr))))
    pks_all <- chromPeaks(x)
    pks_nr <- nrow(pks_all)
    if (length(peakCol) != pks_nr) 
      peakCol <- rep(peakCol[1], pks_nr)
    if (length(peakBg) != pks_nr) 
      peakBg <- rep(peakBg[1], pks_nr)
    if (length(peakPch) != pks_nr) 
      peakPch <- rep(peakPch[1], pks_nr)
    for (i in seq_len(nr)) {
      x_sub <- x[i, , drop = FALSE]
      plot(as(x_sub, ifelse(is(x_sub, "XChromatograms"), 
                            "Chromatograms", "Chromatogram")), col = col, 
           lty = lty, type = type, xlab = xlab, ylab = ylab, 
           main = main, ...)
      idx <- which(pks_all[, "row"] == i)
      if (length(idx) && peakType != "none") {
        pks <- chromPeaks(x_sub)
        .add_chromatogram_peaks(x_sub, pks, col = peakCol[idx], 
                                bg = peakBg[idx], type = peakType, pch = peakPch[idx], 
                                ...)
      }
    }
  }
  .local(x, ...)
}



fig <- fig %>% layout(
  xaxis = list(
    tick0 = min(rt)
  )
)

chrmaken <- {
    .local <- function (object, rt, mz, aggregationFun = "sum", 
                        missing = NA_real_, msLevel = 1L, BPPARAM = bpparam(), 
                        adjustedRtime = hasAdjustedRtime(object), filled = FALSE, 
                        include = c("apex_within", "any", "none")) 
    {
      include <- match.arg(include)
      if (adjustedRtime) 
        adj_rt <- rtime(object, adjusted = TRUE)
      object_od <- as(object, "OnDiskMSnExp")
      object_od <- selectFeatureData(object_od, fcol = c("fileIdx", 
                                                         "spIdx", "seqNum", "acquisitionNum", "msLevel", "polarity", 
                                                         "retentionTime", "precursorScanNum"))
      if (adjustedRtime) {
        object_od@featureData$retentionTime <- adj_rt
      }
      res <- MSnbase::chromatogram(object_od, rt = rt, mz = mz, 
                                   aggregationFun = aggregationFun, missing = missing, 
                                   msLevel = msLevel, BPPARAM = BPPARAM)
      if (!hasChromPeaks(object) | include == "none") 
        return(res)
      lvls <- 1:length(fileNames(object))
      if (missing(rt)) 
        rt <- c(-Inf, Inf)
      if (missing(mz)) 
        mz <- c(-Inf, Inf)
      if (is.matrix(rt) | is.matrix(mz)) {
        if (!is.matrix(rt)) 
          rt <- matrix(rt, ncol = 2)
        if (!is.matrix(mz)) 
          mz <- matrix(mz, ncol = 2)
        if (nrow(rt) == 1) 
          rt <- matrix(rep(rt, nrow(mz)), ncol = 2, byrow = TRUE)
        if (nrow(mz) == 1) 
          mz <- matrix(rep(mz, nrow(rt)), ncol = 2, byrow = TRUE)
        pk_list <- vector("list", nrow(mz))
        pkd_list <- vector("list", nrow(mz))
        for (i in 1:nrow(mz)) {
          pks <- chromPeaks(object, rt = rt[i, ], mz = mz[i, 
          ], type = include)
          pkd <- chromPeakData(object)[rownames(pks), , 
                                       drop = FALSE]
          if (!filled) {
            pks <- pks[!pkd$is_filled, , drop = FALSE]
            pkd <- pkd[!pkd$is_filled, , drop = FALSE]
          }
          smpls <- factor(pks[, "sample"], levels = lvls)
          pk_list[[i]] <- split.data.frame(pks, smpls)
          pkd_list[[i]] <- split.data.frame(pkd, smpls)
        }
        pks <- do.call(rbind, pk_list)
        pks <- pks[seq_along(pks)]
        pkd <- do.call(rbind, pkd_list)
        pkd <- pkd[seq_along(pkd)]
      }
      else {
        pks <- chromPeaks(object, rt = rt, mz = mz, type = include)
        pkd <- chromPeakData(object)[rownames(pks), , drop = FALSE]
        if (!filled) {
          pks <- pks[!pkd$is_filled, , drop = FALSE]
          pkd <- pkd[!pkd$is_filled, , drop = FALSE]
        }
        smpls <- factor(pks[, "sample"], levels = lvls)
        pks <- split.data.frame(pks, smpls)
        pkd <- split.data.frame(pkd, smpls)
      }
      res <- as(res, "XChromatograms")
      res@.Data <- matrix(mapply(unlist(res), pks, pkd, FUN = function(chr, 
                                                                       pk, pd) {
        chr@chromPeaks <- pk
        chr@chromPeakData <- pd
        chr
      }), nrow = nrow(res), dimnames = dimnames(res))
      res@.processHistory <- object@.processHistory
      if (hasFeatures(object)) {
        pks_sub <- chromPeaks(res)
        fts <- lapply(seq_len(nrow(res)), function(r) {
          fdev <- featureDefinitions(object, mz = mz(res)[r, 
          ], rt = rt)
          if (nrow(fdev)) {
            fdev$row <- r
            .subset_features_on_chrom_peaks(fdev, chromPeaks(object), 
                                            pks_sub)
          }
          else DataFrame()
        })
        res@featureDefinitions <- do.call(rbind, fts)
      }
      validObject(res)
      res
    }
    .local(object, ...)
  }
  