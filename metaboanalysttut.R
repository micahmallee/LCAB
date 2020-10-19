library(MetaboAnalystR)
library(googledrive)
data_folder_Sample <- "~/Data_IBD"
data_folder_QC <- "~/QC_IBD"

# Use Google API for data downloading here. 
# Please "install.packages('googledrive')" and "install.packages('httpuv')"first.
library(googledrive);
temp <- tempfile(fileext = ".zip")
# Please authorize your google account to access the data
dl <- drive_download(
  as_id("10DBpPEWy2cZyXvmlLOIfpwqQYwFplKYK"), path = temp, overwrite = TRUE)
# Setting your own date file folder
out <- unzip(temp, exdir = data_folder_QC)
# Date files for parameters optimization are deposited below
out

temp <- tempfile(fileext = ".zip")
dl <- drive_download(
  as_id("1-wlFUkzEwWX1afRWLJlY_KEJs7BfsZim"), path = temp, overwrite = TRUE)
# Setting the date file folder
out <- unzip(temp, exdir = data_folder_Sample)
# Date files for normal processing example are deposited below
out

raw_data <- PerformDataTrimming(data_folder_QC,rt.idx = 0.2)

param_initial <- SetPeakParam(platform = "UPLC-Q/E") 
param_optimized <- PerformParamsOptimization(raw_data, param = param_initial, ncore = 8)

rawData <- ImportRawMSData(data_folder_Sample,plotSettings = SetPlotParam(Plot=F))

mSet <- PerformPeakProfiling(rawData,param_optimized$best_parameters,
                             plotSettings = SetPlotParam(Plot = T))
annParams <- SetAnnotationParam(polarity = "negative", mz_abs_add = 0.005)
annotPeaks <- PerformPeakAnnotation(mSet, annParams)
maPeaks <- FormatPeakList(annotPeaks, annParams, filtIso =F, filtAdducts = FALSE,
                          missPercent = 1)