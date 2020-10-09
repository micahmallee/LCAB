# Readying environment
# Load packages
library(xcms)
library(MSnbase)
library(msdata)
library(magrittr)
library(png)

# Load data
kruid_data <- readMSData(files = "Kruid 131 Zwarte peper 6 191119me_71.mzXML", mode = "onDisk")

head(fData(kruid_data))
rtime(kruid_data)

kruid_chr <- chromatogram(kruid_data)
plot(kruid_chr)

plot(chromatogram(data_cent))
