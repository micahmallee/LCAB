library(MetaboAnalystR)

fls <- dir(path = "mzxml", full.names = TRUE)
pd <- data.frame(file = basename(fls),
                 sample = c("Peper_5", "Peper_6"),
                 group = "Peper")

peperdata <- readMSData(fls, pdata = new("NAnnotatedDataFrame", pd),
                   mode = "onDisk") 


peperdata2 <- ImportRawMSData("/home/admin_bi/Documents/Micah/LCAB/mzxml/")

param_initial <- SetPeakParam(platform = "general")
param_optimized <- PerformParamsOptimization(okee, param = param_initial)

oke <- PerformPeakProfiling(rawData = peperdata, Params = param_initial, plotSettings = SetPlotParam(Plot = TRUE))
