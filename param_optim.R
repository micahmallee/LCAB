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


system.time({
  ## Annotation steps 1 sample
  ### Load data
  spiked_data <- readMSData(files = c('rhino_data/mzxml/Spike/Mara_spike_alkanen.mzXML', 
                                      'rhino_data/mzxml/Spike/Naima_spike_alkanen.mzXML', 
                                      'rhino_data/mzxml/Spike/Vungu_spike_alkanen.mzXML'), mode = 'onDisk')
  spiked_data@phenoData@data[["sample_name"]] <- spiked_data@phenoData@data[["sampleNames"]]
  spiked_data@phenoData@data[["sampleNames"]] <- NULL
  ### Peak detection
  first_eval <- Noise_evaluate(spiked_data)
  ## Define the rt and m/z range of the peak area
  rtr <- c(1840, 1880)
  mzr <- c(334.9, 335.1)
  ## extract the chromatogram
  chr_raw <- chromatogram(spiked_data, mz = mzr, rt = rtr)
  plot(chr_raw)
  
  spiked_data %>%
    filterRt(rt = rtr) %>%
    filterMz(mz = c(72.8, 73.4)) %>%
    plot(type = "XIC")
  # ppm = 70
  # peakwidth = 1,10
  # 
  
  params <- SetPeakParam(ppm = 70, noise = 5000, value_of_prefilter = 5000,
                         prefilter = 2, min_peakwidth = 1, max_peakwidth = 10)
  mSet_spiked_data <- PerformPeakPicking(spiked_data, param = updateRawSpectraParam(params))
  mSet_spiked_data <- PerformPeakAlignment(mSet_spiked_data, param = updateRawSpectraParam(params))
  mSet_spiked_data <- PerformPeakFiling(mSet_spiked_data, param = updateRawSpectraParam(params))
  spiked_data_xcmslist <-  split_mSet(mSet = mSet_spiked_data)
  spiked_data_xcmslist <- annotate_xcmslist(xcmslist = spiked_data_xcmslist, perfwhm = 0.6)
  spiked_data_msp <- lapply(spiked_data_xcmslist, to.msp, file = NULL, settings = NULL, ndigit = 3, minfeat = 10, minintens = 0, intensity = "maxo", secs2mins = F)
  
  ### Get SPLASH hashes and match thirdblocks
  querySPLASH <- get_splashscores(msp_list = spiked_data_msp)
  full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
  full_mona_SPLASHES <- sapply(mona_msp, function(x){
    str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
  })
  query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
  database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
  SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
  similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = spiked_data_msp, database = mona_msp)
  
  ### Get x top matches
  bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.6)
  matchmatrices <- vector(mode = 'list', length = length(bestmatches))
  names(matchmatrices) <- spiked_data_xcms@phenoData[["sample_name"]]
  naampjes <- list()
  bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
  bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
  for (i in seq_along(bestmatches)) {
    matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
    colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
    naampjes[[i]] <- names(bestmatches[[i]])
  }
})