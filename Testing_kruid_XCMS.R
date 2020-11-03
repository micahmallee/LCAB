# Readying environment
# Load packages
library(xcms)
library(MSnbase)
library(msdata)
library(magrittr)
library(png)
library(baseline)

# Load data
xcms_kruid_data <- readMSData(files = "mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML", mode = "onDisk")



head(fData(xcms_kruid_data))
rtime(xcms_kruid_data)

xcms_kruid_chr <- chromatogram(xcms_kruid_data)
plot(xcms_kruid_chr)

filt_kruid_data <- filterRt(kruid_data, c(500, 2500))
plot(chromatogram(filt_kruid_data))
# Findpeaks:
xcms_cwp <- CentWaveParam(peakwidth = c(3, 37.5), ppm = 22.35, snthresh = 5, mzCenterFun = "wMean", prefilter = c(2, 0.02), integrate = 1, mzdiff = 0, fitgauss = F, noise = 0, verboseColumns = F)
xcms_cwp2 <- CentWaveParam(peakwidth = c(3, 37.5), ppm = 22.35, snthresh = 10, mzCenterFun = "wMean", prefilter = c(2, 0.02), integrate = 1, mzdiff = 0, fitgauss = F, noise = 0, verboseColumns = F)
xcms_peaks <- findChromPeaks(xcms_kruid_data, param = xcms_cwp)
xcms_peaks2 <- findChromPeaks(xcms_kruid_data, param = xcms_cwp2)
plot(chromatogram(xcms_peaks2))
# Plot peaks
chromPeaks(xcms_peaks)
plot(chromatogram(xcms_peaks))

# Refine peaks
mpp <- MergeNeighboringPeaksParam(expandRt = 4)
refined_xcms_peaks <- refineChromPeaks(xcms_peaks, param = mpp)
# Show refined peaks
chromPeaks(refined_xcms_peaks)
plot(chromPeaks(refined_xcms_peaks))



# niet cwt maar matched filter
rawxcms <- xcmsRaw(filename = 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML')
xcms_matchedfilter_peaks <- findPeaks.matchedFilter(rawxcms)
xcms_cwt <- findPeaks.centWave(rawxcms)
xcms_foundpeaks <- findPeaks(rawxcms)


# Nu met beide samples:
# Data inladen
pd <- data.frame(file = basename(dir('kruiden/', full.names = T)),
                 sample = c("Kruid_130", "Kruid_131"),
                 group = "Pepers")
xcms_beide <- readMSData(files = dir('kruiden/', full.names = T), mode = 'onDisk', pdata = pd)

# Peak profiling parameters bepalen
xcms_cwp <- CentWaveParam(peakwidth = c(3, 37.5), ppm = 22.35, snthresh = 5, mzCenterFun = "wMean", prefilter = c(2, 0.02), integrate = 1, mzdiff = 0, fitgauss = F, noise = 0, verboseColumns = F)
# Peaks vinden
xcms_peaks_beide <- findChromPeaks(xcms_beide, param = xcms_cwp)

# Peaks visualiseren
chromPeaks(xcms_peaks_beide)
plot(chromatogram(xcms_peaks_beide))

# Refine peaks
mpp <- MergeNeighboringPeaksParam(expandRt = 4)
refined_xcms_peaks_beide <- refineChromPeaks(xcms_peaks_beide, param = mpp)
# Show refined peaks
chromPeaks(refined_xcms_peaks_beide)
plot(chromatogram(refined_xcms_peaks_beide))

pdp <- PeakDensityParam(bw = 2, minFraction = 0.5)
data_cent <- groupChromPeaks(data_cent, pdp)

#' Define settings for the alignment
pgp <- PeakGroupsParam(minFraction = 0.5, span = 0.25)
refined_xcms_peaks_beide <- adjustRtime(refined_xcms_peaks_beide, param = pgp)






oke <- xcmsSet(files = dir('kruiden/', full.names = T), )
xcms::group(oke)
fillPeaks(object = oke)
