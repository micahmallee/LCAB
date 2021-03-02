library(metaMS)
library(MetaboAnalystR)
library(xcms)
library(OrgMassSpecR)

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

## Create pseudospectra
xsetxcms <- xsAnnotate(xs = smSet$xcmsSet)
xsetxcms <- groupFWHM(object = xsetxcms)
xsetxcms2 <- groupCorr(xsetxcms)
### Correct format peakinfo
convertedgroupedxset <- xsetxcms
convertedgroupedxset@groupInfo <- convertedgroupedxset@xcmsSet@peaks
xsetxcmsmsp <- to.msp(object = convertedgroupedxset, file = NULL, settings = metaSetting(TSQXLS.GC, 'DBconstruction'))

DB.treated2 <- treat.DB(mona_msp)


# get mass spectra for each pseudospectrum
#' oke je hebt dus je pseudospectra
#' die geven de peaknummers aan met nagenoeg dezelfde rt
#' om dan het massa spectrum te krijgen moet je op die rt het achterliggende massa spectrum pakken. Deze kan je dan als het goed is gebruiken
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
DB.treated3 <- DB.treated2[1:5]


spectrasimmilarities <- c()
for (i in yee.msp) {
  lapply(1:length(DB.treated2), function(x) {
    spectrasimmilarities <- c(spectrasimmilarities, SpectrumSimilarity(spec.top = yee.msp[[i]][, 1:2], spec.bottom = DB.treated2[[x]]$pspectrum, print.graphic = F, print.alignment = F))
  })
}

result <- matrix(0, length(DB.treated3), length(yee.msp[[1]]))

test <- sapply(1:length(DB.treated3), function(i) {
  result <- SpectrumSimilarity(spec.top = yee.msp[[1]][[1]][, 1:2], spec.bottom = DB.treated3[[i]]$pspectrum, print.alignment = F, output.list = F, print.graphic = F)
  result
})


match.results <- lapply(1:length(yee.msp), function(ii) {
  result <- matrix(0, length(DB.treated3), length(yee.msp[[ii]]))
  for (i in 1:length(yee.msp[[ii]])) {
    result <- lapply(1:length(DB.treated3), function(jj) {
      print(paste0(ii, ' ', i, ' ', jj))
      resultaat <- SpectrumSimilarity(spec.top = yee.msp[[ii]][[i]][, 1:2], spec.bottom = DB.treated3[[jj]]$pspectrum, print.graphic = F, print.alignment = F)
      result[i, jj] <- resultaat
      # print(result)
      result
    })
    print(result)
  }
  result
})


