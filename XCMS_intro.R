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




