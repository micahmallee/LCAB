library(metaMS)
library(MetaboAnalystR)
library(xcms)
library(OrgMassSpecR)
library(splashR)
library(stringr)
library(metaMSdata)
data(FEMsettings)
rm(Orbitrap.RP, Synapt.NP, Synapt.RP)


MSConvert_CMD <- paste0("docker run --rm -v `pwd`:`pwd` -w `pwd` chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ", in_path, " --mzXML")
system(MSConvert_CMD)

test_data <- readMSData('test_data_mzxml/Test_nieuwe_kolom_MCX_SPE_KW1_Spike_18.mzXML', mode = 'onDisk')
trimmed_test_data <- PerformDataTrimming(datapath = c("test_data_mzxml/"), rt.idx = 1, plot = F)
params_opt <- PerformParamsOptimization(raw_data = trimmed_test_data, param = SetPeakParam(platform = 'general', Peak_method = 'centWave'))

p23 <- Noise_evaluate(test_data)
params <- SetPeakParam(ppm = 50, noise = p23$noise, value_of_prefilter = p23$value_of_prefilter, prefilter = p23$prefilter)
smSet_test_data <- PerformPeakPicking(test_data, param = updateRawSpectraParam(params))
smSet_test_data[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
# smSet_test_data <- PerformPeakAlignment(smSet_test_data, param = updateRawSpectraParam(params))
# smSet_test_data <- PerformPeakFiling(smSet_test_data, param = updateRawSpectraParam(params))
oke <- PerformPeakAnnotation(mSet = smSet_test_data, annotaParam = SetAnnotationParam())

oke1 <- new(Class = 'xcmsSet')

ja <- xcms::findChromPeaks(object = test_data, param = CentWaveParam(ppm = 22, peakwidth = c(2, 10), prefilter = c(0, 0) ))
ja1 <- xcms::findChromPeaks(object = test_data, param = CentWaveParam(ppm = 22, peakwidth = c(2, 10), prefilter = c(p23$prefilter, p23$value_of_prefilter) ))
ja2 <- xsAnnotate(ja1)
ja2 <- groupFWHM(ja2)
jamsp <- to.msp(object = ja2,settings = metaSetting(TSQXLS.GC, 'DBconstruction'), file = 'okee.msp')

querysplash <- getsplashscores(list(jamsp))
query_secondblocks <- getsecondblocks(querysplash)
splashmatches <- matchsecondblocks(querysecondblocks = query_secondblocks, databasesecondblocks = mona_secondblocks)
similarity_scores <- similarities(msp_query = list(jamsp), database = mona_msp, SPLASH_hits = splashmatches)


# Compound linking
## Data loading, database, peakprofiling etc
### LCAB:
data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
### Home:
data <- readMSData(files = c('kruiden/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'kruiden/Kruid 126 Zwarte peper 1 191119me_66.mzXML'), mode = 'onDisk')
# Read MoNA database
mona_msp <- read.msp(file = 'MoNA-export-GC-MS_Spectra.msp')
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
# annParams <- SetAnnotationParam(polarity = 'positive')
# # Annotated has pseudospectra (not sure if needed, some steps might be LC-MS specific)
# annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams)


# Split mSet to nested xcmslist for further annotation
annotatedxcmslist <- split_annotate(smSet)

# Convert to MSP format
mspxcmslist <- lapply(X = annotatedxcmslist, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))


kleinedb <- mona_msp[1:100]

full_mona_splashes <- getsplashscores(mona_msp)
mona_secondblocks <- getsecondblocks(full_mona_splashes)

querysplash <- getsplashscores(mspxcmslist)
query_secondblocks <- getsecondblocks(querysplash)

# Match secondblocks
splashmatches <- matchsecondblocks(querysecondblocks = query_secondblocks, databasesecondblocks = mona_secondblocks)

# Get similarity scores
similarity_scores <- similarities(msp_query = mspxcmslist, database = mona_msp, SPLASH_hits = splashmatches)




tophits <- function(similarity_scores, limit, database, splashmatches) {
  topmatches
}



similarities <- function(msp_query, database, SPLASH_hits) {
  #' Calculate Spectrumsimilarity per SPLASH match
  # Loop through hits
  allmatches <- vector(mode = 'list', length = length(SPLASH_hits))
  allmatches <- lapply(seq_along(SPLASH_hits), function(x){
    matches_per_sample <- vector(mode = 'list', length = length(SPLASH_hits[[x]]))
    lapply(seq_along(SPLASH_hits[[x]]), function(y){
      matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(SPLASH_hits[[x]][[y]]))
      matches_per_pseudospectrum <- sapply(SPLASH_hits[[x]][[y]], function(z){
        SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
      })
    })
  })
  # replace NaN values with 0
  allmatches <- lapply(allmatches, function(x){
    x <- lapply(x, function(y){
      y <- sapply(y, function(z){
        z <- ifelse(is.nan(z), 0, z)
      })
    })
  })
  return(allmatches)
}

matchsecondblocks <- function(querysecondblocks, databasesecondblocks) {
  #' This function takes the vectors containing the second SPLASH blocks for all pseudospectra.
  #' Matches second block with second blokcs in database per sample. Returns list with matches per pseudospectra
  matchindexes <- vector(mode = 'list', length = length(querysecondblocks))
  matchindexes <- lapply(querysecondblocks, function(x) {
    matchindexes1 <- vector(mode = 'list', length = length(x))
    matchindexes1 <- sapply(x, function(y){
      grep(pattern = y, x = databasesecondblocks)
    })
  })
}


getsecondblocks <- function(splashscores){
  #' This function takes all second blocks from the SPLASH codes per pseudospectrum.
  if (class(splashscores) == "character") {
    blocks <- vector(mode = 'character', length = length(splashscores))
    blocks <- sapply(splashscores, function(x){
      str_split(string = x, pattern = '-', simplify = F)[[1]][[2]]
    })
  } else {
    blocks <- vector(mode = "list", length = length(splashscores))
    blocks <- lapply(seq_along(splashscores), function(x) {
      blocks_psample <- vector(mode = 'character', length = length(splashscores[[x]]))
      blocks_psample <- sapply(splashscores[[x]], function(y) {
        str_split(string = y, pattern = '-', simplify = F)[[1]][[2]]
      })
    })
  }
  return(blocks)
}

getsplashscores <- function(msp_object) {
  #' This function calculates all SPLASH codes per pseudospectrum
  if (class(msp_object[[1]][[1]]) == "character") {
    splashscoresquery <- vector(mode = "character", length = length(msp_object))
    splashscoresquery <- sapply(seq_along(msp_object), function(y) {
      getSplash(msp_object[[y]]$pspectrum)
    })
  } else {
    splashscoresquery <- vector(mode = "list", length = length(msp_object))
    splashscoresquery <- lapply(seq_along(msp_object), function(x) {
      lapply(1:length(msp_object[[x]]), function(y) {
        getSplash(msp_object[[x]][[y]][, 1:2])
      })
    })
  }
  return(splashscoresquery)
}

split_annotate <- function(mSet) {
  #' This function splits the mSet object based on sample. 
  #' Subsequently this function groups the peaks into pseudospectra.
  f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
  f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
    gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
  })
  splitxcms <- xcms:::split.xcmsSet(smSet$xcmsSet, f = factor(f))
  annotatedxcmslist <- vector(mode = 'list', length = length(f))
  annotatedxcmslist <- lapply(splitxcms, function(x){
    x <- xsAnnotate(x)
    x <- groupFWHM(x)
  })
}