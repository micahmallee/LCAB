library(MetaboAnalystR)
library(googledrive)
data_folder_Sample <- "~/Data_IBD"
data_folder_QC <- "~/QC_IBD"

raw_data1 <- PerformDataTrimming(data_folder_QC,rt.idx = 0.2)

param_initial <- SetPeakParam(platform = "UPLC-Q/E") 
param_optimized1 <- PerformParamsOptimization(raw_data1, param = param_initial, ncore = 1)

rawData <- ImportRawMSData(data_folder_Sample,plotSettings = SetPlotParam(Plot=F))

mSet <- PerformPeakProfiling(rawData,param_optimized,
                             plotSettings = SetPlotParam(Plot = F))
annParams <- SetAnnotationParam(polarity = "negative", mz_abs_add = 0.005)
annotPeaks <- PerformPeakAnnotation(mSet, annParams)
maPeaks <- FormatPeakList(annotPeaks, annParams, filtIso =F, filtAdducts = FALSE,
                          missPercent = 1)

mSet1 <- InitDataObjects("pktable", "stat", FALSE)

mSet1 <- Read.TextData(mSet1, "metaboanalyst_input.csv", "colu", "disc")
mSet2 <- SanityCheckData(mSet)


mSet1 <- ReplaceMin(mSet1)

mSet2 <- Ttests.Anal(mSet, F, 0.25, FALSE, TRUE)

