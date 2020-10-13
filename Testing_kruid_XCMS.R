# Readying environment
# Load packages
library(xcms)
library(MSnbase)
library(msdata)
library(magrittr)
library(png)

# Load data
kruid_data <- readMSData(files = "mzxml/Kruid 131 Zwarte peper 6 191119me_71.mzXML", mode = "onDisk")

head(fData(kruid_data))
rtime(kruid_data)

kruid_chr <- chromatogram(kruid_data)
plot(kruid_chr)

# default centwave parameters:
cwp <- CentWaveParam()

# dry-run peak detection on XIC:
res <- findChromPeaks(kruid_chr, param = cwp)
chromPeaks(res)
plot(res)

cwp <- CentWaveParam(peakwidth = c(2, 10), integrate = 2)
srn_chr <- findChromPeaks(srn_chr, param = cwp)
plot(srn_chr)
chromPeaks(srn_chr)