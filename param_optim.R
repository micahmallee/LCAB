library(MetaboAnalystR)
# library(OptiLCMS)
{
# raw_data1 <- ImportRawMSData(foldername = 'rhino_data/params/', mode = 'inMemory', plotSettings = SetPlotParam(Plot = F))

raw_data1 <- MetaboAnalystR::PerformDataTrimming(datapath = 'rhino_data/params/', plot = F, rt.idx = 1)
param_init <- MetaboAnalystR::SetPeakParam()
param_optim <- MetaboAnalystR::PerformParamsOptimization(raw_data = raw_data1, param = param_init, ncore = 1)


# raw_data2 <- ImportRawMSData(foldername = 'rhino_data/mzxml/', mode = 'onDisk', plotSettings = SetPlotParam(Plot = F))



spiked_data <- readMSData(files = c('rhino_data/mzxml/Spike/Mara_spike_alkanen.mzXML', 
                                    'rhino_data/mzxml/Spike/Naima_spike_alkanen.mzXML', 
                                    'rhino_data/mzxml/Spike/Vungu_spike_alkanen.mzXML'), mode = 'onDisk')

duplo_data_1 <- readMSData(files = c('rhino_data/mzxml/Mara.mzXML', 
                                 'rhino_data/mzxml/Naima.mzXML', 
                                 'rhino_data/mzxml/Vungu.mzXML'), mode = 'onDisk')
duplo_data_2 <- readMSData(files = c('rhino_data/mzxml/Mara_duplo.mzXML', 
                                     'rhino_data/mzxml/Naima_duplo.mzXML', 
                                     'rhino_data/mzxml/Vungu_duplo.mzXML'), mode = 'onDisk')

RColorBrewer::brewer.pal(n = 3, name = 'Pastel1')
tc <- split(tic(duplo_data_1), f = fromFile(duplo_data_1))
boxplot(tc, ylab = "Intensity", main = "Total ion current duplo samples 1")

tc1 <- split(tic(duplo_data_2), f = fromFile(duplo_data_2))
boxplot(tc1, xlab = "Intensity", main = "Total ion current duplo samples 2", names = c('Mara', 'Naima', 'Vungu'), horizontal = T)


oke <- xcmsSet(files = 'rhino_data/mzxml/Spike/Naima_spike_alkanen.mzXML', mslevel = 1)
oke_annotated <- xsAnnotate(xs = oke)
oke_grouped <- groupFWHM(object = oke_annotated, perfwhm = 1)
oke_mps <- lapply(list(oke_grouped), to.msp, file = NULL, settings = NULL, minfeat = 10, minintens = 0.0015, intensity = "maxo", secs2mins = F)
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
}


similarity_scores <- readRDS('shiny/sim_scores')
mSet <- readRDS('shiny/mSet')

ok213e <- lapply(seq_along(bestmatches), function(x){
  top_compounds_per_sample <- sapply(seq_along(bestmatches[[x]]), function(y) {
    bestmatches[[x]][[y]][[1]][1]
  })
})
names(ok213e) <- mSet[["onDiskData"]]@phenoData@data[["sample_name"]]
heatmap_data <- function(sample_compound) {
  jaaa <- sapply(sample_compound, str_c)
  ja <- unique(unlist(jaaa))
  compound_matrix1 <- matrix(nrow = length(sample_compound), ncol = length(ja))
  for (i in seq_along(sample_compound)) {
    for (j in seq_along(sample_compound[[i]])) {
      indexje <- match(sample_compound[[i]][[j]], ja)
      compound_matrix1[i, indexje] <- 1
    }
  }
  compound_matrix1[is.na(compound_matrix1)] <- 0
  colnames(compound_matrix1) <- ja
  rownames(compound_matrix1) <- names(sample_compound)
  return(compound_matrix1)
}

ok213e[[mSet[["onDiskData"]]@phenoData@data[["sample_name"]][[1]]]]

mara_vungu <- ok213e[-c(3:4)]
mara_hm <- heatmap_data(ok213e[1:2])
naima_hm <- heatmap_data(ok213e[3:4])
vungu_hm <- heatmap_data(ok213e[5:6])
all_hm <- heatmap_data(ok213e)

mara_vungu <- heatmap_data(mara_vungu)


naima_mara_hm <- heatmap_data(nenm)

pheatmap::pheatmap(mat = t(mara_hm), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                   color = RColorBrewer::brewer.pal(n = 3, name = "Set1"), main = 'Mara')
pheatmap::pheatmap(mat = t(naima_hm), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                   color = RColorBrewer::brewer.pal(n = 3, name = "Set1"), main = 'Naima')
pheatmap::pheatmap(mat = t(vungu_hm), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                   color = RColorBrewer::brewer.pal(n = 3, name = "Set1"), main = 'Vungu')
pheatmap::pheatmap(mat = t(all_hm), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                   color = RColorBrewer::brewer.pal(n = 3, name = "Set1"), main = "All")

pheatmap::pheatmap(mat = t(naima_mara_hm), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                   color = RColorBrewer::brewer.pal(n = 3, name = "Set1"))

# ax <- list(
#   zeroline = TRUE,
#   showline = TRUE,
#   mirror = "ticks",
#   gridcolor = toRGB("gray50"),
#   gridwidth = 2,
#   zerolinecolor = toRGB("red"),
#   zerolinewidth = 4,
#   linecolor = toRGB("black"),
#   linewidth = 6,
#   titlefont=list(
#     family = "Roman",
#     size = 25
#   )
# )
# 
# ay <- list(
#   zeroline = TRUE,
#   showline = TRUE,
#   mirror = "ticks",
#   gridcolor = toRGB("gray50"),
#   gridwidth = 2,
#   zerolinecolor = toRGB("red"),
#   zerolinewidth = 4,
#   linecolor = toRGB("black"),
#   linewidth = 6
# )
# 
# plot_ly(z = t(mara_hm), type = 'heatmap', x = rownames(mara_hm), 
#         y = colnames(mara_hm), colors = RColorBrewer::brewer.pal(n = 3, name = "Set1"), 
#         hoverinfo = "y") %>% layout(xaxis = ax, yaxis = ay)

system.time({
  ## Annotation steps 1 sample
  ### Load data
  rhino_duplo <- readMSData(files = list.files('rhino_data/duplo_regular', full.names = T), mode = 'onDisk')
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
  
  # params <- SetPeakParam(Peak_method = 'matchedFilter', snthresh = 10, fwhm = 1)
  # mSet_rhino_duplo <- PerformPeakPicking(rhino_duplo, param = updateRawSpectraParam(params))
  # mSet_rhino_duplo <- PerformPeakAlignment(mSet_rhino_duplo, param = updateRawSpectraParam(params))
  # mSet_rhino_duplo <- PerformPeakFiling(mSet_rhino_duplo, param = updateRawSpectraParam(params))
  
  
  
  rhino_duplo_xcmslist <-  split_mSet(mSet = mSet)
  rhino_duplo_xcmslist <- annotate_xcmslist(xcmslist = rhino_duplo_xcmslist, perfwhm = 1)
  rhino_duplo_msp <- lapply(rhino_duplo_xcmslist, to.msp, file = NULL, settings = NULL, minfeat = 10, minintens = 0.0015, intensity = "maxo", secs2mins = F)
  
  ### Get SPLASH hashes and match thirdblocks
  querySPLASH <- get_splashscores(msp_list = rhino_duplo_msp)
  full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
  full_mona_SPLASHES <- sapply(mona_msp, function(x){
    str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
  })
  query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
  database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
  SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
  # similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = rhino_duplo_msp, database = mona_msp)
  
  ### Get x top matches
  bestmatches <- tophits(similarity_scores = similarity_scores, limit = 15, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.65)
  matchmatrices <- vector(mode = 'list', length = length(bestmatches))
  names(matchmatrices) <- names(rhino_duplo_xcmslist)
  naampjes <- list()
  bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
  bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
  for (i in seq_along(bestmatches)) {
    matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
    colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
    naampjes[[i]] <- names(bestmatches[[i]])
  }
})
