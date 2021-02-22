# Testing MetaMS
library(metaMS)
library(metaMSdata)
data(threeStdsDB)
data(FEMsettings)

# Complete metaMS workflow
cdfdir <- system.file('extdata', package = 'metaMSdata')
cdffiles <- list.files(cdfdir, pattern = "_GC_", full.names = T, ignore.case = T)
result <- runGC(files = cdffiles, settings = TSQXLS.GC, DB = DB, nSlaves = 2)

# Whole proces behind runGC:
## Peak picking
GCset <- peakDetection(files = cdffiles, settings = metaSetting(TSQXLS.GC, 'PeakPicking'), convert2list = T, nSlaves = 2)
## Definition of pseudospectra
allSamples <- lapply(GCset, runCAMERA, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA"))
allSamples.msp <- lapply(allSamples, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))
sapply(allSamples.msp, length)
plotPseudoSpectrum(allSamples.msp[[1]][[26]])
## Annotation
DB.treated <- treat.DB(DB)
allSam.matches <- matchSamples2DB(allSamples.msp, DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)
allSam.matches
matchExpSpec(allSamples.msp[[1]][[4]], DB.treated, DB.treated = TRUE, plotIt = TRUE)
## Unknowns
allSamples.msp.scaled <- lapply(allSamples.msp, treat.DB, isMSP = F)
allSamples.matches <- matchSamples2Samples(allSamples.msp.scaled, allSamples.msp, annotations = allSam.matches$annotations, settings = metaSetting(TSQXLS.GC, "betweenSamples"))
names(allSamples.matches)
## Output
features.df <- getFeatureInfo(stdDB = DB, allMatches = allSamples.matches, sampleList = allSamples.msp)
features.df[, c(1:3, ncol(features.df) - 2:0)]
PseudoSpectra <- constructExpPseudoSpectra(allMatches = allSamples.matches, standardsDB = DB)
ann.df <- getAnnotationMat(exp.msp = allSamples.msp, pspectra = PseudoSpectra, allMatches = allSamples.matches)
ann.df
ann.df2 <- sweep(ann.df, 1, sapply(PseudoSpectra, function(x) max(x$pspectrum[, 2])), FUN = "*")
ann.df2
result <- list(PeakTable = cbind(data.frame(features.df), data.frame(round(ann.df2))), PseudoSpectra = PseudoSpectra, settings = TSQXLS.GC, 
     xset = allSamples, annotation = allSam.matches$annotation, samples.msp = allSamples.msp, SessionInfo = sessionInfo())


# Own data:
mmsdb <- read.msp(file = 'MoNA-export-GC-MS_Spectra.msp')
mmssettings <- TSQXLS.GC
## Make data suitable for runCAMERA? or for to.msp
xsetxcms <- xsAnnotate(xs = smSet$xcmsSet)
### Create pseudospectra
groupedxset <- groupFWHM(object = xsetxcms)

# xsetCam <- runCAMERA(smSet$xcmsSet, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA"))

okemsp <- to.msp(groupedxset, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))

okemsp <- result
okematches <- matchSamples2DB(okemsp, DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)



## Definition of pseudospectra
allSamples <- lapply(GCset, runCAMERA, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA"))
allSamples.msp <- lapply(allSamples, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))
sapply(allSamples.msp, length)
plotPseudoSpectrum(allSamples.msp[[1]][[26]])
## Annotation
DB.treated <- treat.DB(DB)
allSam.matches <- matchSamples2DB(allSamples.msp, DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)
allSam.matches
matchExpSpec(allSamples.msp[[1]][[4]], DB.treated, DB.treated = TRUE, plotIt = TRUE)
## Unknowns
allSamples.msp.scaled <- lapply(allSamples.msp, treat.DB, isMSP = F)
allSamples.matches <- matchSamples2Samples(allSamples.msp.scaled, allSamples.msp, annotations = allSam.matches$annotations, settings = metaSetting(TSQXLS.GC, "betweenSamples"))
names(allSamples.matches)
## Output
features.df <- getFeatureInfo(stdDB = DB, allMatches = allSamples.matches, sampleList = allSamples.msp)
features.df[, c(1:3, ncol(features.df) - 2:0)]
PseudoSpectra <- constructExpPseudoSpectra(allMatches = allSamples.matches, standardsDB = DB)
ann.df <- getAnnotationMat(exp.msp = allSamples.msp, pspectra = PseudoSpectra, allMatches = allSamples.matches)
ann.df
ann.df2 <- sweep(ann.df, 1, sapply(PseudoSpectra, function(x) max(x$pspectrum[, 2])), FUN = "*")
ann.df2
result <- list(PeakTable = cbind(data.frame(features.df), data.frame(round(ann.df2))), PseudoSpectra = PseudoSpectra, settings = TSQXLS.GC, 
               xset = allSamples, annotation = allSam.matches$annotation, samples.msp = allSamples.msp, SessionInfo = sessionInfo())



# Convert MAR xcmsSet to xcms xcmsSet, no method for conversion, since slots are the same, loop through.
xset <- smSet$xcmsSet
annxset <- xsAnnotate(xs = xset, sample = c(1:2), polarity = 'positive')
setjesgroeperen <- groupFWHM(object = annxset)

for (i in slotNames(xset)) {
  slot(lekkersetjehoor, i) <- slot(xset, i)
}

lijstje <- list()
lijstje[[1]] <- xs


object = xsetCam
xsetCam@groupInfo <- smSet$xcmsSet@peaks


function (object, file = NULL, settings = NULL, ndigit = 0, minfeat, 
          minintens, intensity = c("maxo", "into"), secs2mins = TRUE) 
{
  if (!is.null(settings)) {
    intensity <- settings$intensityMeasure
    minfeat <- settings$minfeat
    minintens <- settings$minintens
  }
  else {
    intensity <- match.arg(intensity)
  }
  if (class(object) == "xsAnnotate") {
    allpks <- object@groupInfo
    minI <- minintens * max(allpks[, intensity])
    tooSmall <- which(allpks[, intensity] < minI)
    pspectra <- lapply(object@pspectra, function(x) x[!x %in% 
                                                        tooSmall])
  }
  else {
    minI <- minintens * max(object[, intensity])
    allpks <- object[object[, intensity] >= minI, ]
    pspectra <- split(1:nrow(allpks), allpks[, "rt"])
  }
  npeaks <- sapply(pspectra, length)
  pspectra <- pspectra[npeaks >= minfeat]
  if (!is.null(file)) {
    if (length(pspectra) > 0) {
      for (i in 1:length(pspectra)) {
        ofile <- paste(file, "_", ceiling(i/1000), ".txt", 
                       sep = "")
        newfile <- (i%%1000) == 1
        idx <- pspectra[[i]]
        pks <- allpks[idx, , drop = FALSE]
        pks <- pks[order(pks[, "mz"]), , drop = FALSE]
        pks[, intensity] <- 1000 * pks[, intensity]/max(pks[, 
                                                            intensity])
        cat("Name: grp ", i, " (rt: ", mean(pks[, "rt"]), 
            ")", sep = "", file = ofile, append = !newfile)
        cat("\nNum Peaks:", nrow(pks), file = ofile, 
            append = TRUE)
        for (ii in 1:nrow(pks)) cat("\n(", round(pks[ii, 
                                                     "mz"], ndigit), "\t", round(pks[ii, intensity], 
                                                                                 ndigit), ")", file = ofile, append = TRUE, 
                                    sep = "")
        cat("\n\n", file = ofile, append = TRUE)
      }
    }
  }
  result <- metaMS:::removeDoubleMasses(lapply(pspectra, function(x) cbind(mz = round(allpks[x, 
                                                                                    "mz"], digits = ndigit), allpks[x, c(intensity, "rt")])))
  if (secs2mins) {
    invisible(lapply(result, function(x) cbind(x[, c("mz", 
                                                     intensity)], rt = x[, "rt"]/60)))
  }
  else {
    invisible(result)
  }
}
