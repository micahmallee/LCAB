library(MetaboAnalystR)
library(OptiLCMS)

# raw_data1 <- ImportRawMSData(foldername = 'rhino_data/params/', mode = 'inMemory', plotSettings = SetPlotParam(Plot = F))
# raw_data2 <- ImportRawMSData(foldername = 'rhino_data/mzxml/', mode = 'onDisk', plotSettings = SetPlotParam(Plot = F))

# Sys.setenv('R_REMOTES_NO_ERRORS_FROM_WARNINGS' = 'true')

raw_data_trimmed <- OptiLCMS::PerformROIExtraction(datapath = 'rhino_mzxml/Regular/', plot = F, rmConts = F)
param_initial <- OptiLCMS::SetPeakParam(platform = 'general', Peak_method = 'centWave', snthresh = 10, max_peakwidth = 16, min_peakwidth = 1)
param_optimized <- OptiLCMS::PerformParamsOptimization(mSet = raw_data_trimmed, param = param_initial, ncore = 4)
saveRDS(param_optimized, 'optimized_params')
# param_optimized <- PerformParamsOptimization(raw_data, param = param_initial, ncore = 1)
# 


# test_data_trimmed <- OptiLCMS::PerformROIExtraction(datapath = 'test_data_mzxml/converted/', rmConts = F, rt.idx = 1)
# test_data_inmem <- readMSData('test_data_mzxml/raw/spike18.mzXML', mode = 'inMemory', msLevel. = 1)
# param_optimized <-OptiLCMS::PerformParamsOptimization(test_data_trimmed, param = OptiLCMS::SetPeakParam(noise = 10, min_peakwidth = 2, max_peakwidth = 10), ncore = 1)
# 
# oke_dan <- PerformPeakPicking(object = oke, param = SetPeakParam(Peak_method = 'centWave', snthresh = 10, min_peakwidth = 2, max_peakwidth = 10))
# 
# oke <- OptiLCMS::ImportRawMSData(path = c('test_data_mzxml/raw/spike18.mzXML'), mode = 'onDisk')
