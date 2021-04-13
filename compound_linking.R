library(metaMS)
library(MetaboAnalystR)
library(xcms)
library(OrgMassSpecR)
library(splashR)
library(stringr)
# data(FEMsettings)
# rm(Orbitrap.RP, Synapt.NP, Synapt.RP)

in_path <- 'c://Users/Micah/Documents/LCAB/test_data_mzxml/raw/spike18.raw'
MSConvert_CMD <- paste0("docker run --rm -v `pwd`:`pwd` -w `pwd` chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ", in_path, " --mzXML")
system(MSConvert_CMD)

tweesamples <- readMSData(files = list.files('kruiden/', full.names = T), mode = 'onDisk')


plot_chrom_tic_bpc <- function(raw_data) {
  plotData <- data.frame(scantime = rtime(raw_data), tic = tic(raw_data), bpc = bpi(raw_data))
  p <- plot_ly(source = "p") %>% 
    add_trace(data = plotData, x = ~scantime, y = ~tic, type = "scatter", mode = "lines", line = list(color = "rgba(0, 0, 0,0.7)", width = 0.8), 
              text = ~paste(scantime, "s"), name = "<b>Total ion chromatogram</b>") %>% 
    add_trace(data = plotData, x = ~scantime, y = ~bpc, type = "scatter", mode = "lines", line = list(color = "rgba(0, 215, 167,1)", width = 1.1), 
              text = ~paste(round(bpc, 4), "m/z"), name = "<b>Base peak chromatogram</b>") %>% 
    layout(legend = list(x = 0.7, y = 0.99), 
           xaxis = list(title = "Scan time (s)", range = c(0, max(plotData$scantime)), showspikes = TRUE, spikemode = "toaxis+across", spikesnap = "data", 
                        showline = FALSE, zeroline = FALSE, spikedash = "solid", showgrid = TRUE), 
           yaxis = list(title = "Counts", showgrid = FALSE, showticklabels = TRUE, zeroline = FALSE, showline = FALSE), hovermode = "x", showlegend = TRUE) %>% 
    event_register("plotly_click")
  return(p)
}


oke <- plot_chrom_tic_bpc(tweesamples)



# Preload MoNA_DB and SPLASH hashes
mona_msp <- readRDS(file = 'shiny/data/mona_msp')
mona_splashes <- readRDS(file = 'shiny/data/mona_splashes')
system.time({
  ## Annotation steps 1 sample
  ### Load data
  test_data <- readMSData('test_data_mzxml/Test_nieuwe_kolom_MCX_SPE_KW1_Spike_18.mzXML', mode = 'onDisk')
  # mona_msp <- read.msp('MoNA-export-GC-MS_Spectra.msp')
  # trimmed_test_data <- PerformDataTrimming(datapath = c("test_data_mzxml/"), rt.idx = 1, plot = F)
  # params_opt <- PerformParamsOptimization(raw_data = trimmed_test_data, param = SetPeakParam(platform = 'general', Peak_method = 'centWave'))
  ### Peak detection
  p23 <- Noise_evaluate(test_data)
  params <- SetPeakParam(ppm = 22, noise = p23$noise, value_of_prefilter = p23$value_of_prefilter, prefilter = p23$prefilter, min_peakwidth = 1, max_peakwidth = 15, snthresh = 10, mzdiff = 0.1)
  smSet_test_data <- PerformPeakPicking(test_data, param = updateRawSpectraParam(params))
  smSet_test_data[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]]
  smSet_test_data[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
  
  ### Get pseudospectra and convert to msp format
  test_data_xcms <- mSet2xcmsSet(smSet_test_data)
  test_data_xsannotate <- xsAnnotate(test_data_xcms)
  test_data_xsannotate <- groupFWHM(test_data_xsannotate, perfwhm = 1)
  test_data_msp <- to.msp(object = test_data_xsannotate, file = NULL, settings = NULL, ndigit = 3, minfeat = 5, minintens = 0, intensity = "maxo", secs2mins = F)
  
  ### Get SPLASH hashes and match thirdblocks
  querySPLASH <- get_splashscores(msp_list = list(test_data_msp))
  full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
  full_mona_SPLASHES <- sapply(mona_msp, function(x){
    str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
  })
  query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
  database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
  SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
  similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = test_data_msp, database = mona_msp)
  
  ### Get x top matches
  bestmatches <- tophits(similarity_scores = similarity_scores, limit = 15, database = mona_msp, splashmatches = SPLASH_matches)
  names(x = bestmatches) <- c(seq_along(bestmatches))
  bestmatches <- Filter(function(k) length(k) > 0, bestmatches)
  matchmatrix <- t(data.table::rbindlist(bestmatches))
})



## Multiple samples
msamples <- readMSData(c('mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML', 'mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
p24 <- Noise_evaluate(msamples)
params2 <- SetPeakParam(ppm = p24$ppm, noise = p24$noise, value_of_prefilter = p24$value_of_prefilter, prefilter = p24$prefilter, min_peakwidth = 2, max_peakwidth = 8, snthresh = 100, mzdiff = 0.1)
smSet_msamples <- PerformPeakPicking(msamples, param = updateRawSpectraParam(params2))
smSet_msamples[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet_msamples[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet_msamples[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet_msamples <- PerformPeakAlignment(smSet_msamples, param = updateRawSpectraParam(params2))
smSet_msamples <- PerformPeakFiling(smSet_msamples, param = updateRawSpectraParam(params2))
xcmslist <-  split_annotate(mSet = smSet_msamples)
msplist <- lapply(xcmslist, to.msp, file = NULL, settings = NULL, ndigit = 3, minfeat = 5, minintens = 0, intensity = "maxo", secs2mins = F)

## Compare multiple sample approach to single sample approach
blabla <- readMSData('kruiden/Kruid 126 Zwarte peper 1 191119me_66.mzXML', mode = 'onDisk')
blabla1 <- PerformPeakPicking(blabla, param = updateRawSpectraParam(params2))
blabla1[["onDiskData"]]@phenoData@data[["sample_name"]] <- blabla1[["onDiskData"]]@phenoData@data[["sampleNames"]]
blabla1[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
blablaxcms <- mSet2xcmsSet(blabla1)
blablaxsannotate <- xsAnnotate(blablaxcms)
blablaxsannotate <- groupFWHM(blablaxsannotate, perfwhm = 3)
blablamsp <- to.msp(object = blablaxsannotate, file = NULL, settings = NULL, ndigit = 3, minfeat = 5, minintens = 0, intensity = "maxo", secs2mins = F)


for (i in 1:length(bestmatches)) {
  try(expr = {
    if (bestmatches[[i]][[1]][["Score"]] > 80) {
      print(paste0(names(bestmatches)[i], '  Best match: ', bestmatches[[i]][[1]][[1]][["Name"]], ' Score: ', bestmatches[[i]][[1]][["Score"]]))
    }
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
      index1 <- grep(pattern = y, x = unlist(similarity_scores[[1]][[x]]))
      if (!is.na(splashmatches[[1]][[x]][index1]) & as.numeric(y) > 0.8) {
        paste0('Compound: ', as.character(database[[splashmatches[[1]][[x]][index1[1]]]]["Name"]), '   Score: ', round((y * 100), 2))
        # matchandscore[["Match"]] <- as.character(database[[splashmatches[[1]][[x]][index1[1]]]]["Name"])
        # matchandscore[["Score"]] <- round((y * 100), 2)
      } else {
        # matchandscore <- NULL
        NULL
      }
      # matchandscore
    })
  })
}

similarities_thirdblocks <- function(nine_matches, msp_query, database) {
  #' Calculate Spectrumsimilarity per SPLASH match
  # Loop through hits
  if (length(nine_matches) == 1) {
    msp_query <- list(msp_query)
  }
  if (.Platform$OS.type == 'windows') {
    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    clusterExport(cl=cl, 'mona_msp')
    allmatches <- vector(mode = 'list', length = length(nine_matches))
    allmatches <- lapply(seq_along(nine_matches), function(x) {
      matches_per_sample <- vector(mode = 'list', length = length(nine_matches[[x]]))
      matches_per_sample <- parLapply(cl, seq_along(nine_matches[[x]]), function(y){
        matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(nine_matches[[x]][[y]]))
        matches_per_pseudospectrum <- sapply(nine_matches[[x]][[y]], function(z){
          OrgMassSpecR::SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
        })
      })
    })
    stopCluster(cl)
  } else {
    num_cores <- detectCores()
    allmatches <- vector(mode = 'list', length = length(nine_matches))
    allmatches <- lapply(seq_along(nine_matches), function(x) {
      matches_per_sample <- vector(mode = 'list', length = length(nine_matches[[x]]))
      matches_per_sample <- mclapply(mc.cores = num_cores, seq_along(nine_matches[[x]]), function(y){
        matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(nine_matches[[x]][[y]]))
        matches_per_pseudospectrum <- sapply(nine_matches[[x]][[y]], function(z){
          OrgMassSpecR::SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
        })
      })
    })
  }
  return(allmatches)
}

match_nines <- function(query_blocks, database_blocks) {
  #' This function takes the vectors containing the third SPLASH blocks for all pseudospectra.
  #' Matches the position of the character 9 in the the third blocks in database per sample. 
  #' Returns list with matches per pseudospectra.
  db_nines <- str_locate(database_blocks, '9')
  db_nines <- db_nines[, 1]
  query_nines <- str_locate(query_blocks, '9')
  query_nines <- query_nines[, 1]
  indexes_p_sample <- vector(mode = 'list', length = length(query_nines))
  indexes_p_sample <- lapply(seq_along(query_nines), function(y) {
    indexes <- sapply(seq_along(db_nines), function(z) ifelse(test = query_nines[[y]] == db_nines[[z]], yes = z, no = 0))
    indexes <- indexes[indexes != 0]
  })
  return(indexes_p_sample)
}


get_blocks <- function(splashscores, blocknr = 3){
  #' This function takes all n blocks from the SPLASH codes per pseudospectrum.
  if (class(splashscores) == "character") {
    blocks <- str_split(string = splashscores, pattern = '-')
    blocks <- sapply(blocks, function(x) x[[blocknr]])
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


get_splashscores <- function(msp_list) {
  #' This function calculates all SPLASH codes per pseudospectrum
  totalsplashscores <- vector(mode = 'list', length = length(msp_list))
  totalsplashscores <- lapply(msp_list, function(x) {
    splashscoresquery <- vector(mode = "character", length = length(x))
    splashscoresquery <- sapply(x, function(y) getSplash(y[, 1:2]))
  })
  return(totalsplashscores)
}


mSet2xcmsSet <-  function(mSet) {
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
  xs@mslevel <- 1
  xs@scanrange <- range(scanIndex(mSet[[format]]))
  if (any(mSet$msFeatureData$chromPeakData$is_filled)) {
    fld <- which(mSet$msFeatureData$chromPeakData$is_filled)
    xs@filled <- as.integer(fld)
  }
  return(xs)
}

split_annotate <- function(mSet) {
  #' This function splits the mSet object based on sample. 
  #' Subsequently this function groups the peaks into pseudospectra.
  f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
  f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
    gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
  })
  splitxcms <- xcms:::split.xcmsSet(mSet$xcmsSet, f = factor(f))
  annotatedxcmslist <- vector(mode = 'list', length = length(f))
  annotatedxcmslist <- lapply(splitxcms, function(x){
    x <- xsAnnotate(x)
    x <- groupFWHM(x, perfwhm = 3)
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
