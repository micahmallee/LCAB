library(MetaboAnalystR)
library(xcms)


fls <- dir(path = "mzxml", full.names = TRUE)
pd <- data.frame(file = basename(fls),
                 sample = c("Peper_5", "Peper_6"),
                 group = "Peper")

peperdata <- readMSData(fls, pdata = new("NAnnotatedDataFrame", pd),
                   mode = "onDisk") 




peperdata2 <- ImportRawMSData("mzxml/")



oke <- PerformPeakProfiling(rawData = peperdata, Params = param_initial, plotSettings = SetPlotParam(Plot = TRUE))

kruiden <- PerformDataTrimming("mzxml/", rt.idx = 1)
param_initial <- SetPeakParam(platform = "general")
param_optimized <- PerformParamsOptimization(raw_data = kruiden, param = param_initial)
mSet <- PerformPeakProfiling(rawData = kruiden, Params = param_optimized)
annParams <- SetAnnotationParam(polarity = "positive")
annotPeaks <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams)
maPeaks <- FormatPeakList(annotPeaks = annotPeaks, annParams = annParams, filtIso = F, filtAdducts = F, missPercent = 1)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "metaboanalyst_input.csv", "colu", "disc")




raw_kruid_6 <- PerformDataTrimming("mzxml/", rt.idx = 1)
param_initial <- SetPeakParam(platform = "general")
param_optimized <- PerformParamsOptimization(raw_kruid_6, param = param_initial)

mSet <- PerformPeakProfiling(raw_kruid_6, Params = param_initial)
