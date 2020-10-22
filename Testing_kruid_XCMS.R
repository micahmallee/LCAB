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

filt_kruid_data <- filterRt(kruid_data, c(500, 2500))
plot(chromatogram(filt_kruid_data))
# even een peak eruit halen
<<<<<<< HEAD
cwp <- CentWaveParam(peakwidth = c(10, 30))
okeee <- findChromPeaks(filt_kruid_data, param = cwp)
=======
>>>>>>> 33a24d16e580fe338d6631274b55e8c25ecacce5


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