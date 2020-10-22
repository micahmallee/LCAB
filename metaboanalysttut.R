library(MetaboAnalystR)
library(googledrive)
data_folder_Sample <- "~/Data_IBD"
data_folder_QC <- "~/QC_IBD"

raw_data <- PerformDataTrimming(data_folder_QC,rt.idx = 0.2)

param_initial1 <- SetPeakParam(platform = "UPLC-Q/E") 
param_optimized <- PerformParamsOptimization(raw_data, param = param_initial1, ncore = 8)

rawData <- ImportRawMSData(data_folder_Sample,plotSettings = SetPlotParam(Plot=F))

mSet <- PerformPeakProfiling(rawData,param_optimized$best_parameters,
                             plotSettings = SetPlotParam(Plot = T))
annParams <- SetAnnotationParam(polarity = "negative", mz_abs_add = 0.005)
annotPeaks <- PerformPeakAnnotation(mSet, annParams)
maPeaks <- FormatPeakList(annotPeaks, annParams, filtIso =F, filtAdducts = FALSE,
                          missPercent = 1)