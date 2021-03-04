library(metaMS)
library(MetaboAnalystR)
library(xcms)
library(OrgMassSpecR)
library(splashR)


# Compound linking
## Data loading, database, peakprofiling etc
### LCAB:
data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
### Home:
data <- readMSData(files = c('kruiden/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'kruiden/Kruid 126 Zwarte peper 1 191119me_66.mzXML'), mode = 'onDisk')
param_optimized <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0.01,
                           snthresh = 10, bw = 2, ppm = 20, min_peakwidth = 5, max_peakwidth = 30, 
                           noise = 100, prefilter = 3, value_of_prefilter = 100, minFraction = 0.5, 
                           minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                           family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                           mzCenterFun = "wMean")

# Create OnDiskData and msFeatureData
smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
# Adds to onDiskData and creates and fills FeatureGroupTable
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized))
# Creates and fills xcmsSet
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized))

# Annotation parameters and create annotated smSet
annParams <- SetAnnotationParam(polarity = 'positive')
# Annotated has pseudospectra (not sure if needed, some steps might be LC-MS specific)
annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams)

# Read MoNA database
mona_msp <- read.msp(file = 'MoNA-export-GC-MS_Spectra.msp')


mona_splashes <- getsplashscores(mona_msp)
kleinedb <- mona_msp[1:5]

oke <- similarities(msp_query = yee.msp, database = mona_msp)


getsplashscores <- function(msp_object) {
  splashscoresquery <- vector(mode = "list", length = length(msp_object))
  splashscoresquery <- lapply(1: length(msp_object), function(x) {
    lapply(1:length(msp_object[[x]]), function(y) {
      getSplash(msp_object[[x]][[y]][, 1:2])
    })
  })
  return(splashscoresquery)
}


result <- matrix(0, length(database), length(msp_query[[1]]))
result <- lapply(1:length(msp_query[[1]]), function(y) {
  scores <- vector(mode = "numeric", length = length(database))
  scores <- sapply(database, function(z) {
    SpectrumSimilarity(spec.top = msp_query[[1]][[1]][, 1:2], spec.bottom = z$pspectrum, print.graphic = F, print.alignment = F)
  })
})


result <- matrix(0, length(database), length(msp_query[[1]]))
scores <- vector(mode = "numeric", length = length(msp_query[[1]]))

dbhit <- sapply(1:length(database), function(z) {
  scores <- vector(mode = "numeric", length = length(database))
  scores[z] <- SpectrumSimilarity(spec.top = msp_query[[1]][[1]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.graphic = F, print.alignment = F)
  scores
})
result[z, y] <- SpectrumSimilarity(spec.top = msp_query[[1]][[1]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.graphic = F, print.alignment = F)



similarities <- function(msp_query, database) {
  scores <- vector(mode = "numeric", length = length(database))
  spectrasimmilarities <- vector(mode = "list", length = length(msp_query))
  spectrasimmilarities <- lapply(1:length(msp_query), function(x) {
    result <- matrix(0, length(database), length(msp_query[[x]]))
    result <- lapply(1:length(msp_query[[x]]), function(y) {
      scores <- sapply(database, function(z) {
        SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = z$pspectrum, print.graphic = F, print.alignment = F)
      })
    })
    result
  })
}





f <- c()
for (i in smSet$xcmsSet@phenoData$sample_name) {f <- (c(f, gsub(x = i, " ", ".", fixed = T)))}
yee <- xcms:::split.xcmsSet(smSet$xcmsSet, f = factor(f))
yee2 <- c()
for (i in yee) {
  oke <- xsAnnotate(xs = i)
  oke <- groupFWHM(object = oke)
  yee2 <- c(yee2, oke)
}
yee.msp <- lapply(yee2, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))


