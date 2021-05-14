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

rhino <- readMSData(files = list.files('rhino_data/params/', full.names = T), mode = 'inMemory', msLevel. = 1)
rhino_trimmed <- PerformDataTrimming( datapath = 'rhino_data/params', plot = SetPlotParam(Plot = F), write = T, rt.idx = 1)
raw <- ImportRawMSData(foldername = 'rhino_data/params', mode = 'inMemory', plotSettings = SetPlotParam(Plot = F))

# Plot all peaks per sample
p <- plot_chrom_tic_bpc(smSet_msamples$onDiskData, tic_visibility = 'legendonly')
for (i in seq_along(xcmslist)) {
  p <- p %>% add_trace(data = data.frame(xcmslist[[i]]@groupInfo), 
                       x = ~rt, y = ~maxo, type = "scatter", mode = "markers", text = ~mz, 
                       name = paste(xcmslist[[i]]@xcmsSet@phenoData$sample_name))
}

# Plot all pseudospectra per sample
plotData_single_sample <- plotdata_pseudospectra(list(test_data_msp))
plotData_single_sample1 <- data.table::rbindlist(plotData, use.names = T)

plotData_msamples <- plotdata_pseudospectra(msplist)
plotData_msamples1 <- data.table::rbindlist(plotData_msamples, use.names = T)

p <- plot_chrom_tic_bpc(smSet_msamples$onDiskData, tic_visibility = 'legendonly')
p <- plotdata_pseudospectra_traces(plotData_msamples, xcmslist = xcmslist, p)

# Plot all found compounds
plotData_compounds <- get_compound_plotData(bestmatches = bestmatches1, plotData = plotData_msamples)
p <- plot_chrom_tic_bpc(smSet_msamples$onDiskData, tic_visibility = 'legendonly')
p <- plotdata_compounds_traces(plotData = plotData_compounds, p = p, msplist = msplist)

# single sample
plotData_compounds <- get_compound_plotData(bestmatches = bestmatches, plotData = plotdata_pseudospectra(test_data_msp))
p <- plot_chrom_tic_bpc(smSet_test_data$onDiskData, tic_visibility = 'legendonly', source = 'p2')
p <- plotdata_compounds_traces(plotData = plotData_compounds, p = p, msplist = test_data_msp)

plotdata_compounds_traces <- function(plotData, p, msplist) {
  for (i in seq_along(plotData)) {
    p <- p %>% add_trace(data = plotData[[i]], x = ~rt, y = ~intense, type = "scatter", mode = "markers", text = ~Compound, name = ~paste0("Compounds: ", names(msplist)[i]))
  }
  return(p)
}

plotdata_pseudospectra_traces <- function(plotData, xcmslist, p){
  for (i in seq_along(plotData)) {
    p <- p %>% add_trace(data = plotData[[i]], x = ~rt, y = ~intense, type = "scatter", mode = "markers", text = ~pspectra_id, name = ~paste0('Pseudospectra: ', xcmslist[[i]]@xcmsSet@phenoData$sample_name))
  }
  return(p)
}


get_plotData_pseudospectra <- function(msps){
  plotData <- vector(mode = 'list', length = length(msps))
  for (i in seq_along(msps)){
    pseudospectrum <- lapply(seq_along(msps[[i]]), function(x){
      rt <- round(median(msps[[i]][[x]][, 3]), 4)
      intense <- max(msps[[i]][[x]][, 2])
      list(rt = rt, intense = intense, pspectra_id = x, sample_nr = i)
    })
    plotData[[i]] <- pseudospectrum
  }
  plotData <- lapply(plotData, data.table::rbindlist, use.names = T)
  return(plotData)
}


oke <- vector(mode = 'list', length = length(xcmslist))
for (i in seq_along(xcmslist)) {
  oke[[i]] <- as.data.frame(xcmslist[[i]]@peaks)
  }
oke1 <- rbind(oke[[1]], oke[[2]])
oke1$rt <- round(oke1$rt, 4)


get_plotData_compounds <- function(bestmatches, plotData) {
  compound_plotData <- vector('list', length = length(bestmatches))
  compound_plotData <- lapply(seq_along(bestmatches), function(x) {
    pspectra_nr <- str_extract(string = names(bestmatches[[x]]), "[0-9]")
    compound_plotData_per_sample <- plotData[[x]][as.integer(pspectra_nr),]
    top_compounds_per_sample <- vector(mode = 'character', length = length(bestmatches[[x]]))
    top_compounds_per_sample <- sapply(seq_along(bestmatches[[x]]), function(y) {
      bestmatches[[x]][[y]][[1]][1]
    })
    compound_plotData_per_sample <- cbind(compound_plotData_per_sample, Compound = top_compounds_per_sample)
    compound_plotData_per_sample
  })
}


plot_align_spectra <- function(topSpectrum, botSpectrum) {
  top_spectrum <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
  top_spectrum$normalized <- round((top_spectrum$intensity/max(top_spectrum$intensity)) * 100)
  fig <- plot_ly() %>% add_bars(
    y = -1*abs(top_spectrum$normalized),
    x = top_spectrum$mz,
    name = 'Pseudospectrum sample'
  )
  fig <- fig %>% add_bars(
    x = spec.bottom[,1],
    y = spec.bottom[,2],
    name = 'Pseudospectrum database'
  )
  return(fig)
}


# Preload MoNA_DB and SPLASH hashes
mona_msp <- readRDS(file = 'shiny/data/mona_msp')
mona_splashes <- readRDS(file = 'shiny/data/mona_splashes')
system.time({
  ## Annotation steps 1 sample
  ### Load data
  # test_data <- readMSData('test_data_mzxml/Test_nieuwe_kolom_MCX_SPE_KW1_Spike_18.mzXML', mode = 'onDisk')
  test_data <- readMSData('test_data_mzxml/converted/spike18.mzXML', mode = 'onDisk')
  test_data@phenoData@data[["sample_name"]] <- test_data@phenoData@data[["sampleNames"]]
  test_data@phenoData@data[["sampleNames"]] <- NULL
  # mona_msp <- read.msp('MoNA-export-GC-MS_Spectra.msp')
  # trimmed_test_data <- PerformDataTrimming(datapath = c("test_data_mzxml/"), rt.idx = 1, plot = F)
  # params_opt <- PerformParamsOptimization(raw_data = trimmed_test_data, param = SetPeakParam(platform = 'general', Peak_method = 'centWave'))
  ### Peak detection
  # p23 <- Noise_evaluate(test_data)
  params <- SetPeakParam(ppm = 22, noise = 10000, value_of_prefilter = 0.01, prefilter = 2, min_peakwidth = 1, max_peakwidth = 10, snthresh = 100)
  smSet_test_data <- PerformPeakPicking(test_data, param = updateRawSpectraParam(params))
  test_data_xcms <- mSet2xcmsSet(smSet_test_data)
  
  ### Get pseudospectra and convert to msp format
  test_data_xsannotate <- xsAnnotate(test_data_xcms)
  test_data_xsannotate <- groupFWHM(test_data_xsannotate, perfwhm = 0.6)
  test_data_msp <- to.msp(object = test_data_xsannotate, file = NULL, settings = NULL, ndigit = 0, minfeat = 3, minintens = 0, intensity = "maxo", secs2mins = F)
  
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
  bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.6)
  matchmatrices <- vector(mode = 'list', length = length(bestmatches))
  names(matchmatrices) <- test_data_xcms@phenoData[["sample_name"]]
  naampjes <- list()
  bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
  bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
  for (i in seq_along(bestmatches)) {
    matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
    colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
    naampjes[[i]] <- names(bestmatches[[i]])
  }
})



## Multiple samples
msamples <- readMSData(c('mzxml/KRUID_126/Kruid 126 Zwarte peper 1 191119me_66.mzXML', 'mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
msamples@phenoData@data[["sample_name"]] <- msamples@phenoData@data[["sampleNames"]]
msamples@phenoData@data[["sampleNames"]] <- NULL
p24 <- Noise_evaluate(msamples)
params2 <- SetPeakParam(ppm = p24$ppm, noise = 10000, value_of_prefilter = 10000, prefilter = p24$prefilter, min_peakwidth = 2, max_peakwidth = 8, snthresh = 1000, mzdiff = 0.1)
smSet_msamples <- PerformPeakPicking(msamples, param = updateRawSpectraParam(params2))
smSet_msamples <- PerformPeakAlignment(smSet_msamples, param = updateRawSpectraParam(params2))
smSet_msamples <- PerformPeakFiling(smSet_msamples, param = updateRawSpectraParam(params2))
xcmslist <-  split_mSet(mSet = smSet_msamples)
xcmslist <- annotate_xcmslist(xcmslist = xcmslist, perfwhm = 0.6)
msplist <- lapply(xcmslist, to.msp, file = NULL, settings = NULL, ndigit = 3, minfeat = 10, minintens = 0, intensity = "maxo", secs2mins = F)
querySPLASH1 <- get_splashscores(msp_list = msplist)
full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
full_mona_SPLASHES <- sapply(mona_msp, function(x){
  str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
})
query_thirdblocks1 <- lapply(querySPLASH1, get_blocks, blocknr = 3)
database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
SPLASH_matches1 <- lapply(query_thirdblocks1, match_nines, database_blocks = database_thirdblocks)
similarity_scores1 <- similarities_thirdblocks(nine_matches = SPLASH_matches1, msp_query = msplist, database = mona_msp)
for(i in seq_along(similarity_scores1)) {names(similarity_scores1[[i]]) <- sprintf("Pseudospectrum_%d", seq.int(similarity_scores1[[i]]))}
bestmatches1 <- tophits(similarity_scores = similarity_scores1, limit = 5, database = mona_msp, splashmatches = SPLASH_matches1, score_cutoff = 0.8)
matchmatrices <- vector(mode = 'list', length = length(bestmatches1))
names(matchmatrices) <- names(xcmslist)
naampjes <- list()
for (i in seq_along(bestmatches1)) {
  matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches1[[i]]))
  colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches1[[i]]), times = 2))
  naampjes[[i]] <- names(bestmatches1[[i]])
}

plotData_peaks <- list(xcmslist[[1]]@peaks)
for (i in 2:length(xcmslist)) {
  plotData_peaks[[i]] <- xcmslist[[i]]@peaks
}

plotData_peaks %>% filter(sample == 2) -> ja

p <- plot_chrom_tic_bpc(smSet_msamples$onDiskData, tic_visibility = 'legendonly', source = 'p')
plotData_peaks <- list(xcmslist[[1]]@peaks)
for (i in 2:length(xcmslist)) {
  plotData_peaks[[i]] <- xcmslist[[i]]@peaks
}
samplename <- sapply(xcmslist, function(x) x@phenoData[["sample_name"]])

for (i in seq_along(plotData_peaks)) {
  p <- p %>% add_trace(data = data.frame(plotData_peaks[[i]]), 
                       x = ~rt, y = ~maxo, type = "scatter", mode = "markers", text = ~mz, 
                       name = paste(samplename[i], ' total peaks'))
}











{
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
}

{
plot_chrom_tic_bpc <- function(raw_data, tic_visibility = NULL, source = NULL) {
  if (is.null(tic_visibility)) {
    hovermode <- "x"
  } else {
    hovermode <- "closest"
  }
  samplenames <- raw_data@phenoData@data[["sample_name"]]
  fdataPerSample <- split.data.frame(raw_data@featureData@data, raw_data@featureData@data$fileIdx)
  cols <- RColorBrewer::brewer.pal(n = length(fdataPerSample) * 2, 'Paired')
  p <- plot_ly(source = source)
  for (i in seq_along(fdataPerSample)) {
    colindex <- if(i%%2 == 0) c(-i+1, -i) else c(i, i+1)
    plotData <- data.frame(scantime = fdataPerSample[[i]]$retentionTime, tic = fdataPerSample[[i]]$totIonCurrent, bpc = fdataPerSample[[i]]$basePeakIntensity)
    p <- p %>% 
      add_trace(data = plotData, x = ~scantime, y = ~tic, type = "scatter", mode = "lines", line = list(color = cols[colindex[1]], width = 0.8), 
                text = ~paste(scantime, "s"), visible = tic_visibility, name = paste("<b>Total ion chromatogram</b> ", samplenames[i])) %>% 
      add_trace(data = plotData, x = ~scantime, y = ~bpc, type = "scatter", mode = "lines", line = list(color = cols[colindex[2]], width = 1.1), 
                text = ~paste(round(bpc, 4), "m/z"), name = paste("<b>Base peak chromatogram</b> ", samplenames[i]))
  }
  p <- p %>% 
    layout(legend = list(x = 0.7, y = 0.99), 
           xaxis = list(title = "Scan time (s)", range = c(0, max(plotData$scantime)), showspikes = TRUE, spikemode = "toaxis+across", spikesnap = "data", 
                        showline = FALSE, zeroline = FALSE, spikedash = "solid", showgrid = TRUE), 
           yaxis = list(title = "Counts", showgrid = FALSE, showticklabels = TRUE, zeroline = FALSE, showline = FALSE), hovermode = "closest", showlegend = TRUE) %>% 
    event_register("plotly_click")
  return(p)
}
  
  
tophits <- function(similarity_scores, limit = 5, database, splashmatches, score_cutoff = 0.8) {
  totalmatches <- vector('list', length(similarity_scores))
  totalmatches <- lapply(seq_along(similarity_scores), function(z){
    indexes_per_sample <- vector(mode = 'list', length = length(similarity_scores[[z]]))
    indexes_per_sample <- lapply(seq_along(similarity_scores[[z]]), function(x){
      top5_scores <- sort(unlist(similarity_scores[[z]][[x]]), decreasing = T)[1:limit]
      top5_matches <- vector(mode = 'list', length = length(top5_scores))
      top5_matches <- lapply(top5_scores, function(y){
        index1 <- grep(pattern = y, x = unlist(similarity_scores[[z]][[x]]))
        if (!is.na(splashmatches[[z]][[x]][index1]) & as.numeric(y) > score_cutoff) {
          c(as.character(database[[splashmatches[[z]][[x]][index1[1]]]]["Name"]), round((y * 100), 2), database[[splashmatches[[z]][[x]][index1[1]]]]["pspectrum"])
        } else {
          NULL
        }
      })
    })
    names(indexes_per_sample) <- c(sprintf("Pseudospectrum_%d", seq.int(indexes_per_sample)))
    indexes_per_sample <- Filter(Negate(function(x) is.null(unlist(x))), indexes_per_sample)
  })
  return(totalmatches)
}


similarities_thirdblocks <- function(nine_matches, msp_query, database) {
  #' Calculate SpectrumSimilarity per SPLASH match
  # Loop through hits
  if (.Platform$OS.type == 'windows') {
    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    clusterExport(cl=cl, c('database', 'msp_query'), envir = environment())
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

split_mSet <- function(mSet) {
  #' This function splits the mSet object based on sample. 
  f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
  f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
    gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
  })
  splitxcms <- xcms:::split.xcmsSet(mSet$xcmsSet, f = factor(f))
  for(i in seq_along(splitxcms)) {splitxcms[[i]]@peaks[, 11] <- i}
  return(splitxcms)
}


annotate_xcmslist <- function(xcmslist, perfwhm = 0.6){
  #' This function groups the peaks into pseudospectra per sample.
  annotatedxcmslist <- vector(mode = 'list', length = length(xcmslist))
  annotatedxcmslist <- lapply(xcmslist, function(x){
    x <- xsAnnotate(x)
    x <- groupFWHM(x, perfwhm = perfwhm)
  })
  return(annotatedxcmslist)
}
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
