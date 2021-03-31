library(metaMS)
library(MetaboAnalystR)
library(xcms)
library(OrgMassSpecR)
library(splashR)
library(stringr)
data(FEMsettings)
rm(Orbitrap.RP, Synapt.NP, Synapt.RP)


MSConvert_CMD <- paste0("docker run --rm -v `pwd`:`pwd` -w `pwd` chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ", in_path, " --mzXML")
system(MSConvert_CMD)

test_data <- readMSData('test_data_mzxml/Test_nieuwe_kolom_MCX_SPE_KW1_Spike_18.mzXML', mode = 'onDisk')
trimmed_test_data <- PerformDataTrimming(datapath = c("test_data_mzxml/"), rt.idx = 1, plot = F)
params_opt <- PerformParamsOptimization(raw_data = trimmed_test_data, param = SetPeakParam(platform = 'general', Peak_method = 'centWave'))

params_pre <- SetPeakParam(platform = 'general', Peak_method = 'centWave')


p23 <- Noise_evaluate(test_data)
params <- SetPeakParam(ppm = p23$ppm, noise = p23$noise, value_of_prefilter = p23$value_of_prefilter, prefilter = p23$prefilter, min_peakwidth = 2, max_peakwidth = 10)
smSet_test_data <- PerformPeakPicking(test_data, param = updateRawSpectraParam(params))
smSet_test_data[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL

test_data_xcms <- mSet2xcmsSet(smSet_test_data)
test_data_xsannotate <- xsAnnotate(test_data_xcms)
test_data_xsannotate <- groupFWHM(test_data_xsannotate)
test_data_msp <- to.msp(object = test_data_xsannotate, file = NULL, settings = NULL, ndigit = 3, minfeat = 1, minintens = 0, intensity = "maxo", secs2mins = F)
test_data_msp <- to.msp(object = ja2,settings = metaSetting(TSQXLS.GC, 'DBconstruction'), file = NULL, intensity = "maxo")

querysplash <- getsplashscores(list(test_data_msp))
query_secondblocks <- getsecondblocks(querysplash)
splashmatches <- matchsecondblocks(querysecondblocks = query_secondblocks, databasesecondblocks = mona_secondblocks)
similarity_scores <- similarities(msp_query = list(test_data_msp), database = mona_msp, SPLASH_hits = splashmatches)
bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = splashmatches)

for (i in 1:length(bestmatches2)) {
  try(expr = {
    print(paste0(i, '  Best match: ', bestmatches2[[i]][[1]][[1]][["Name"]], ' Score: ', bestmatches2[[i]][[1]][["Score"]]))
  }, silent = T)
}


#' Spiked data:
#' L-Valine, TMS derivative       andere splash
#' Glycine, 3TMS derivative       half, wel .86 similarity, via MONA website wel!
#' Serine, 3TMS derivative        Nee
#' L-Threonine, 3TMS derivative   Nee
#' L-Hydroxyproline, (E)-, 3TMS derivative    Ja (36)
#' DL-Phenylalanine, TMS derivative         ja 
#' Phenylalanine, 2TMS derivative   
#' Octopamine, 3TMS derivative      Denk het wel, zit niet in MoNA
#' L-Lysine, 4TMS derivative        OUI MONSIEUR
#' Octopamine, 4TMS derivative      Denk het wel, zit niet in MoNA
#' Palmitic Acid, TMS derivative    Ja zit erin (andere splash code?)
#' Serotonin, 4TMS derivative       ja 5-Hydoxytryptamine

oke <- vector(length(test1[[1]]), mode = 'numeric')
oke2 <- vector(length(test1[[1]]), mode = 'numeric')
for (i in 1:length(test1[[1]])) {
  woord <- str_split(test1[[1]][i], pattern = ',')
  woord[[1]][2]
  oke[[i]] <- as.numeric(woord[[1]][1])
  oke2[[i]] <- as.numeric(woord[[1]][2])
}

oke3 <- rbind(oke, oke2)
oke3 <- t(oke3)
oke3 <- data.frame(oke3)

L_valine_TMS <- oke3
Glycine_TMS <- oke3
Serine_TMS <- oke3
L_Threonine_TMS <- oke3


DL_Phenylalanine_TMS <- oke3
Phenylalanine_TMS <- oke3
Palmatic_acid_TMS <- oke3



L_valine_results <- compare_single_compound(L_valine_TMS, test_data_msp)
L_valine_results2 <- compare_single_compound(l_valine_msp[[1]]$pspectrum, test_data_msp)
Glycine_results <- compare_single_compound(Glycine_TMS, test_data_msp)
Serine_results <- compare_single_compound(Serine_TMS, test_data_msp)
Serine_results2 <- compare_single_compound(serine_msp[[1]]$pspectrum, test_data_msp)
Serine_results3 <- compare_single_compound(serine_msp2[[1]]$pspectrum, test_data_msp)
L_Threonine_results <- compare_single_compound(compound = L_Threonine_TMS, test_data_msp)

l_hydroxyproline_results <- compare_single_compound(compound = l_hydroxyproline[[1]]$pspectrum, database = test_data_msp)

l_lysine_results <- compare_single_compound(l_lysine_msp[[1]]$pspectrum, test_data_msp)

DL_Phenylalanine_results <- compare_single_compound(DL_Phenylalanine_TMS, test_data_msp)
Phenylalanine_results <- compare_single_compound(compound = Phenylalanine_TMS, database = test_data_msp)
Phenylalanine_results2 <- compare_single_compound(compound = Phenylalanine_TMS_2[[1]]$pspectrum, database = test_data_msp)

Palmatic_acid_results <- compare_single_compound(Palmatic_acid_TMS, database = test_data_msp)

silat_results <- compare_single_compound(compound = silat[[1]]$pspectrum, database = jamsp)



compare_single_compound <- function(compound, database) {
  # compoundspectrum <- new('Spectrum2', mz = compound[, 1], intensity = compound[, 2])
  results <- vector(mode = 'list', length = length(database))
  results <- lapply(database, function(x) {
    SpectrumSimilarity(spec.top =  x[, 1:2], spec.bottom = compound, print.alignment = F, print.graphic = F)
    # checkspectrum <- new('Spectrum2', mz = x[, 1], intensity = x[, 2])
    # compareSpectra(checkspectrum, compoundspectrum)
  })
  return(results)
}


oke <- c(rbind(jamsp[[1]][, 1], jamsp[[1]][, 2]))


oke <- c(rbind(test_data_msp[[14]][, 1], test_data_msp[[14]][, 2]))




function (spec.top, spec.bottom, t = 0.25, b = 10, top.label = NULL, 
          bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, 
          print.alignment = FALSE, print.graphic = TRUE, output.list = FALSE) 
{
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 
                                                                 2])
  top_tmp$normalized <- round((top_tmp$intensity/max(top_tmp$intensity)) * 
                                100)
  top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= 
                      xlim[2])
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)
  top <- subset(top_plot, top_plot$intensity >= b)
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 
                                                                          2])
  bottom_tmp$normalized <- round((bottom_tmp$intensity/max(bottom_tmp$intensity)) * 
                                   100)
  bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & 
                         bottom_tmp$mz <= xlim[2])
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)
  bottom <- subset(bottom_plot, bottom_plot$intensity >= b)
  for (i in 1:nrow(bottom)) top[, 1][bottom[, 1][i] >= top[, 
                                                           1] - t & bottom[, 1][i] <= top[, 1] + t] <- bottom[, 
                                                                                                              1][i]
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  if (length(unique(alignment[, 1])) != length(alignment[, 
                                                         1])) 
    warning("the m/z tolerance is set too high")
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  if (print.alignment == TRUE) {
    print(alignment)
  }
  if (x.threshold < 0) 
    stop("x.threshold argument must be zero or a positive number")
  alignment <- alignment[alignment[, 1] >= x.threshold, ]
  u <- alignment[, 2]
  v <- alignment[, 3]
  similarity_score <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))
  if (print.graphic == TRUE) {
    plot.new()
    plot.window(xlim = xlim, ylim = c(-125, 125))
    ticks <- c(-100, -50, 0, 50, 100)
    for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 
                                               2), c(0, top_plot$intensity[i]), col = "blue")
    for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 
                                                  2), c(0, -bottom_plot$intensity[i]), col = "red")
    axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], 
         ylab = "intensity")
    axis(1, pos = -125)
    lines(xlim, c(0, 0))
    rect(xlim[1], -125, xlim[2], 125)
    mtext("m/z", side = 1, line = 2)
    mtext("intensity (%)", side = 2, line = 2)
    plot.window(xlim = c(0, 20), ylim = c(-10, 10))
    text(10, 9, top.label)
    text(10, -9, bottom.label)
  }
  if (output.list == TRUE) {
    headTailPlot <- function() {
      pushViewport(plotViewport(c(5, 5, 2, 2)))
      pushViewport(dataViewport(xscale = xlim, yscale = c(-125, 
                                                          125)))
      grid.rect()
      tmp <- pretty(xlim)
      xlabels <- tmp[tmp >= xlim[1] & tmp <= xlim[2]]
      grid.xaxis(at = xlabels)
      grid.yaxis(at = c(-100, -50, 0, 50, 100))
      grid.segments(top_plot$mz, top_plot$intensity, top_plot$mz, 
                    rep(0, length(top_plot$intensity)), default.units = "native", 
                    gp = gpar(lwd = 0.75, col = "blue"))
      grid.segments(bottom_plot$mz, -bottom_plot$intensity, 
                    bottom_plot$mz, rep(0, length(bottom_plot$intensity)), 
                    default.units = "native", gp = gpar(lwd = 0.75, 
                                                        col = "red"))
      grid.abline(intercept = 0, slope = 0)
      grid.text("intensity (%)", x = unit(-3.5, "lines"), 
                rot = 90)
      grid.text("m/z", y = unit(-3.5, "lines"))
      popViewport(1)
      pushViewport(dataViewport(xscale = c(0, 20), yscale = c(-10, 
                                                              10)))
      grid.text(top.label, unit(10, "native"), unit(9, 
                                                    "native"))
      grid.text(bottom.label, unit(10, "native"), unit(-9, 
                                                       "native"))
      popViewport(2)
    }
    p <- grid.grabExpr(headTailPlot())
  }
  if (output.list == FALSE) {
    return(similarity_score)
  }
  if (output.list == TRUE) {
    return(list(similarity.score = similarity_score, alignment = alignment, 
                plot = p))
  }
}





# Compound linking
## Data loading, database, peakprofiling etc
### LCAB:
data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
### Laptop:
data <- readMSData(files = c('mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', 'mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML'), mode = 'onDisk')
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

full_mona_splashes <- vector(mode = 'character', length = length(mona_msp))
full_mona_splashes <- sapply(mona_msp, function(x){
  str_split(x$Comments, pattern = '\"')[[1]][22]
})



mona_secondblocks <- getsecondblocks(full_mona_splashes)

querysplash <- getsplashscores(mspxcmslist)
query_secondblocks <- getsecondblocks(querysplash)

# Match secondblocks
splashmatches <- matchsecondblocks(querysecondblocks = query_secondblocks, databasesecondblocks = mona_secondblocks)

# Get similarity scores
similarity_scores <- similarities(msp_query = mspxcmslist, database = mona_msp, SPLASH_hits = splashmatches)

# Get top x matches
bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = splashmatches)


tophits <- function(similarity_scores, limit = 5, database, splashmatches) {
  indexes_per_sample <- vector(mode = 'list', length = length(similarity_scores))
  indexes_per_sample[[1]] <- lapply(seq_along(similarity_scores[[1]]), function(x){
    top5_scores <- sort(unlist(similarity_scores[[1]][[x]]), decreasing = T)[1:limit]
    top5_matches <- vector(mode = 'list', length = length(top5_scores))
    top5_matches <- lapply(top5_scores, function(y){
      matchandscore <- list(Match = '', Score = 0)
      index1 <- grep(pattern = y, x = unlist(similarity_scores[[1]][[x]]))
      if (!is.na(splashmatches[[1]][[x]][index1])) {
        matchandscore[["Match"]] <- database[[splashmatches[[1]][[x]][index1[1]]]]
        matchandscore[["Score"]] <- round((y * 100), 2)
      } else {
        matchandscore <- NULL
      }
      matchandscore
    })
  })
}


################### FIXEN matches per sample wordt niet gebruikt hahahaah



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

match_secondblocks <- function(querysecondblocks, databasesecondblocks) {
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

similarities_thirdblocks <- function(msp_query, database, nine_matches) {
  #' Calculate Spectrumsimilarity per SPLASH match
  # Loop through hits
  allmatches <- vector(mode = 'list', length = length(nine_matches))
  allmatches <- lapply(seq_along(nine_matches), function(x){
    matches_per_sample <- vector(mode = 'list', length = length(nine_matches[[x]]))
    matches_per_sample <- lapply(seq_along(nine_matches[[x]]), function(y){
      matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(nine_matches[[x]][[y]]))
      matches_per_pseudospectrum <- sapply(nine_matches[[x]][[y]], function(z){
        SpectrumSimilarity(spec.top = msp_query[[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
      })
    })
  })
  # replace NaN values with 0
  # allmatches <- lapply(allmatches, function(x){
  #   x <- lapply(x, function(y){
  #     y <- sapply(y, function(z){
  #       z <- ifelse(is.nan(z), 0, z)
  #     })
  #   })
  # })
  return(allmatches)
}

match_nines <- function(query_blocks, database_blocks) {
  #' This function takes the vectors containing the third SPLASH blocks for all pseudospectra.
  #' Matches the position of the character 9 in the the third blocks in database per sample. 
  #' Returns list with matches per pseudospectra.
  db_nines <- str_locate(database_blocks, '9')
  query_nines <- str_locate(query_blocks[[1]], '9')
  
  total_indexes <- vector(mode = 'list', length = length(query_nines[, 1]))
  
  total_indexes <- lapply(seq_along(query_nines[, 1]), function(x) {
    indexes <- vector(mode = 'numeric', length = length(db_nines[, 1]))
    
    indexes <- sapply(seq_along(db_nines[, 1]), function(y) ifelse(test = query_nines[,1][[x]] == db_nines[, 1][[y]], yes = y, no = 0))
    
    indexes <- indexes[indexes != 0]
  })
  return(total_indexes)
}

unlist(lapply(twee, function(x) which(een %in% x)))

get_blocks <- function(splashscores, blocknr = 3){
  #' This function takes all third blocks from the SPLASH codes per pseudospectrum.
  if (class(splashscores) == "character") {
    blocks <- vector(mode = 'character', length = length(splashscores))
    blocks <- sapply(splashscores, function(x){
      str_split(string = x, pattern = '-', simplify = F)[[1]][[blocknr]]
    })
  } else {
    blocks <- vector(mode = "list", length = length(splashscores))
    blocks <- lapply(seq_along(splashscores), function(x) {
      blocks_psample <- vector(mode = 'character', length = length(splashscores[[x]]))
      blocks_psample <- sapply(splashscores[[x]], function(y) {
        str_split(string = y, pattern = '-', simplify = F)[[1]][[blocknr]]
      })
    })
  }
  return(blocks)
}


get_splashscores <- function(msp_object) {
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

mSet2xcmsSet <-  function(mSet = smSet_test_data) {
  xs <- new("xcmsSet")
  xs@peaks <- mSet[["msFeatureData"]][["chromPeaks"]]
  rts <- list()
  if (class(mSet[[2]]) == "OnDiskMSnExp") {
    format <- "onDiskData"
  }
  else {
    format <- "inMemoryData"
  }
  rts$raw <- rtime(mSet[[format]])
  rts$corrected <- mSet[["msFeatureData"]][["adjustedRT"]]
  xs@rt <- rts
  xs@phenoData <- pData(mSet[[format]])
  xs@filepaths <- fileNames(mSet[[format]])
  profMethod <- NA
  profStep <- NA
  profParam <- list()
  xs@mslevel <- 1
  xs@scanrange <- range(scanIndex(mSet[[format]]))
  if (any(mSet$msFeatureData$chromPeakData$is_filled)) {
    fld <- which(mSet$msFeatureData$chromPeakData$is_filled)
    xs@filled <- as.integer(fld)
  }
  return(xs)
}