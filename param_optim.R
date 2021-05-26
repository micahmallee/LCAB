library(MetaboAnalystR)
# library(OptiLCMS)

# raw_data1 <- ImportRawMSData(foldername = 'rhino_data/params/', mode = 'inMemory', plotSettings = SetPlotParam(Plot = F))

raw_data1 <- MetaboAnalystR::PerformDataTrimming(datapath = 'rhino_data/params/', plot = F, rt.idx = 1)
param_init <- MetaboAnalystR::SetPeakParam()
param_optim <- MetaboAnalystR::PerformParamsOptimization(raw_data = raw_data1, param = param_init, ncore = 1)


# raw_data2 <- ImportRawMSData(foldername = 'rhino_data/mzxml/', mode = 'onDisk', plotSettings = SetPlotParam(Plot = F))



spiked_data <- readMSData(files = c('rhino_data/mzxml/Spike/Mara_spike_alkanen.mzXML', 
                                    'rhino_data/mzxml/Spike/Naima_spike_alkanen.mzXML', 
                                    'rhino_data/mzxml/Spike/Vungu_spike_alkanen.mzXML'), mode = 'onDisk')
tc <- split(tic(spiked_data), f = fromFile(spiked_data))
boxplot(tc, ylab = "intensity", main = "Total ion current")

oke <- xcmsSet(files = 'rhino_data/mzxml/Spike/Naima_spike_alkanen.mzXML', mslevel = 1)
oke_annotated <- xsAnnotate(xs = oke)
oke_grouped <- groupFWHM(object = oke_annotated, perfwhm = 1)
oke_mps <- lapply(list(oke_grouped), to.msp, file = NULL, settings = NULL, ndigit = 3, minfeat = 10, minintens = 0, intensity = "maxo", secs2mins = F)
querySPLASH <- get_splashscores(msp_list = oke_mps)
full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
full_mona_SPLASHES <- sapply(mona_msp, function(x){
  str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
})
query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = oke_mps, database = mona_msp)

### Get x top matches
bestmatches <- tophits(similarity_scores = similarity_scores, limit = 15, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.6)
matchmatrices <- vector(mode = 'list', length = length(bestmatches))
names(matchmatrices) <- spiked_data_xcmslist@phenoData[["sample_name"]]
naampjes <- list()
bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
for (i in seq_along(bestmatches)) {
  matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
  colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
  naampjes[[i]] <- names(bestmatches[[i]])
}


compoundmatrix <- matrix()

ok213e <- lapply(seq_along(bestmatches), function(x){
  top_compounds_per_sample <- sapply(seq_along(bestmatches[[x]]), function(y) {
    bestmatches[[x]][[y]][[1]][1]
  })
})
jaaa <- sapply(ok213e, str_c)
ja <- unique(unlist(jaaa))

compoundmatrix <- data.frame(row.names = names(bestmatches))

for (i in seq_along(ok213e)) {
  for (j in seq_along(ok213e[[i]])) {
    indexje <- match(ok213e[[i]][[j]], ja)
    compoundmatrix[i, indexje] <- 1
  }
}
compoundmatrix[is.na(compoundmatrix)] <- 0
colnames(compoundmatrix) <- ja

hittemapje <- ggplot(data = compoundmatrix)


system.time({
  ## Annotation steps 1 sample
  ### Load data
  rhino_duplo <- readMSData(files = list.files('rhino_data/duplo_regular/', full.names = T), mode = 'onDisk')
  rhino_duplo@phenoData@data[["sample_name"]] <- rhino_duplo@phenoData@data[["sampleNames"]]
  rhino_duplo@phenoData@data[["sampleNames"]] <- NULL
  ### Peak detection
  # first_eval <- Noise_evaluate(rhino_duplo)
  ## Define the rt and m/z range of the peak area
  # rtr <- c(1840, 1880)
  # mzr <- c(334.9, 335.1)
  ## extract the chromatogram
  # chr_raw <- chromatogram(rhino_duplo, mz = mzr, rt = rtr)
  # plot(chr_raw)
  # 
  # rhino_duplo %>%
  #   filterRt(rt = rtr) %>%
  #   filterMz(mz = c(72.8, 73.4)) %>%
  #   plot(type = "XIC")
  # ppm = 70
  # peakwidth = 1,10
  # 
  
  params <- SetPeakParam(Peak_method = 'matchedFilter', snthresh = 10, fwhm = 1)
  mSet_rhino_duplo <- PerformPeakPicking(rhino_duplo, param = updateRawSpectraParam(params))
  mSet_rhino_duplo <- PerformPeakAlignment(mSet_rhino_duplo, param = updateRawSpectraParam(params))
  mSet_rhino_duplo <- PerformPeakFiling(mSet_rhino_duplo, param = updateRawSpectraParam(params))
  
  rhino_duplo_xcmslist <-  split_mSet(mSet = mSet_rhino_duplo)
  rhino_duplo_xcmslist <- annotate_xcmslist(xcmslist = rhino_duplo_xcmslist, perfwhm = 0.6)
  rhino_duplo_msp <- lapply(rhino_duplo_xcmslist, to.msp, file = NULL, settings = NULL, minfeat = 10, minintens = 0.001, intensity = "maxo", secs2mins = F)
  
  ### Get SPLASH hashes and match thirdblocks
  querySPLASH <- get_splashscores(msp_list = rhino_duplo_msp)
  full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
  full_mona_SPLASHES <- sapply(mona_msp, function(x){
    str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
  })
  query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
  database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
  SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
  similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = rhino_duplo_msp, database = mona_msp)
  
  ### Get x top matches
  bestmatches <- tophits(similarity_scores = similarity_scores, limit = 10, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.6)
  matchmatrices <- vector(mode = 'list', length = length(bestmatches))
  names(matchmatrices) <- rhino_duplo_xcmslist@phenoData[["sample_name"]]
  naampjes <- list()
  bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
  bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
  for (i in seq_along(bestmatches)) {
    matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
    colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
    naampjes[[i]] <- names(bestmatches[[i]])
  }
})
