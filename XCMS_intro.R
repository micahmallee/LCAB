# https://jorainer.github.io/metabolomics2018/xcms-preprocessing.html#34_Preprocessing_of_LC-MS_data
# Install the Bioconductor package manager
install.packages("BiocManager")

# Install the required packages
BiocManager::install(c("xcms",
                       "MSnbase",
                       "msdata",
                       "magrittr",
                       "png"))

# Load packages
library(xcms)
library(MSnbase)
library(msdata)
library(magrittr)
library(png)



## Define the file names.
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)

## Define a data.frame with additional information on the files.
pd <- data.frame(file = basename(fls),
                 injection_idx = c(1, 19),
                 sample = c("POOL_1", "POOL_2"),
                 group = "POOL")

data <- readMSData(fls, pdata = new("NAnnotatedDataFrame", pd),
                   mode = "onDisk") 

#' Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(2)))
} else {
  register(bpstart(SnowParam(2)))
}

pData(data)
data$injection_idx

#' Access spectrum header information
head(fData(data))

# Get retentiontimes
head(rtime(data))

# Get retentiontime split by file
rts <- split(rtime(data), fromFile(data))

# rts is list met length 2. aantal spectra per file: 931 931
lengths(rts)

# Alle spectra tussen 180 en 181 seconden worden opgehaald:
sps <- data %>% filterRt(rt = c(180, 181)) %>% spectra
# de hoeveelheid is de lengte: 6
length(sps)
# met fromFile kan je het herkomstbestand terugvinden:
sapply(sps, fromFile)

# Plot the last mass spectrum in sps:
plot(sps[[6]])

# chromatographic data
chr <- chromatogram(data)
chr
plot(chr)

ints <- intensity(chr[1,1])
head(ints)

pData(chr)

#' Je kan de plot voor serine in dit geval opvragen door te filteren op rt en mz
#' en dit vervolgens te plotten:
data %>%
  filterRt(rt = c(175, 189)) %>%
  filterMz(mz = c(106.02, 106.07)) %>%
  chromatogram(aggregationFun = "max") %>%
  plot()

#' Profile vs centroid mode: profile bevat alle discrete metingen van de continue
#' metingen, centroid reduceert deze metingen tot een betrouwbaar middelpunt.
#' Dit scheelt data

## Centroiden van data:

data %>%
  filterRt(rt = c(175, 189)) %>%
  filterMz(mz = c(106.02, 106.07)) %>%
  plot(type = "XIC") 

#' Smooth the signal, then do a simple peak picking.
data_cent <- data %>%
  smooth(method = "SavitzkyGolay", halfWindowSize = 6) %>%
  pickPeaks()

#' Plot the centroided data for Serine
data_cent %>%
  filterRt(rt = c(175, 189)) %>%
  filterMz(mz = c(106.02, 106.07)) %>%
  plot(type = "XIC") 

#' OMdat het onDiks is gelezen, wordt de data niet aangepast, om het persistent
#' te maken moet de data worden weggeschreven en daarna eventueel weer worden gelezen

#' Write the centroided data to files with the same names in the current
#' directory
fls_new <- basename(fileNames(data))
writeMSData(data_cent, file = fls_new)

#' Read the centroided data.
data_cent <- readMSData(fls_new, pdata = new("NAnnotatedDataFrame", pd),
                        mode = "onDisk") 



srn_chr <- chromatogram(data_cent, rt = c(164, 200),
                        mz = c(106.03, 106.06),
                        aggregationFun = "max")

par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(srn_chr)


# default centwave parameters:
cwp <- CentWaveParam()

# dry-run peak detection on XIC:
res <- findChromPeaks(srn_chr, param = cwp)
chromPeaks(res)

cwp <- CentWaveParam(peakwidth = c(2, 10), integrate = 2)
srn_chr <- findChromPeaks(srn_chr, param = cwp)
plot(srn_chr)
chromPeaks(srn_chr)

#' Restrict the data to signal from Serine
srn <- data_cent %>%
  filterRt(rt = c(179, 186)) %>%
  filterMz(mz = c(106.04, 106.06))

#' Plot the data
plot(srn, type = "XIC") 

#' Extract intensities and split them by file. This will return
#' a list of lists.
ints_per_file <- split(intensity(srn), fromFile(srn))

#' For each file, sum up the intensities.
ints_sum <- lapply(ints_per_file, function(x) sum(unlist(x)))
ints_sum

#' Extract the Serine data for one file as a data.frame
srn_df <- as(filterFile(srn, file = which.max(ints_sum)), "data.frame")

#' The difference between m/z values from consecutive scans expressed
#' in ppm
diff(srn_df$mz) * 1e6 / mean(srn_df$mz) 

# Perform peak detection
cwp <- CentWaveParam(peakwidth = c(2, 10), ppm = 30, integrate = 2)
data_cent <- findChromPeaks(data_cent, param = cwp)

#' Access the peak detection results from a specific m/z - rt area
chromPeaks(data_cent, mz = c(106, 107), rt = c(150, 190)) 

eic_serine <- chromatogram(data_cent, mz = c(106.04, 106.06),
                           rt = c(179, 186))
chromPeaks(eic_serine)

plot(eic_serine)


srn <- data_cent %>%
  filterRt(rt = c(175, 188)) %>%
  filterMz(mz = c(106.04, 106.06))

plot(srn, type = "XIC")

mpp <- MergeNeighboringPeaksParam(expandRt = 4)
data_cent_pp <- refineChromPeaks(data_cent, param = mpp)

chromPeakData(data_cent_pp)
mzr <- c(124.084, 124.088)
rtr <- c(150, 170)
chr_1 <- chromatogram(filterFile(data_cent, 2), mz = mzr, rt = rtr)
chr_2 <- chromatogram(filterFile(data_cent_pp, 2), mz = mzr, rt = rtr)
par(mfrow = c(1, 2))
plot(chr_1)
plot(chr_2)

data_cent <- data_cent_pp

par(mfrow = c(1, 2))
plotChromPeaks(data_cent, 1)
plotChromPeaks(data_cent, 2) 

