# Testing MetaMS
library(metaMS)
library(metaMSdata)
data(threeStdsDB)
data(FEMsettings)

cdfdir <- system.file('extdata', package = 'metaMSdata')
cdffiles <- list.files(cdfdir, pattern = "_GC_", full.names = T, ignore.case = T)
result <- runGC(files = cdffiles, settings = TSQXLS.GC, DB = DB, nSlaves = 2)

testresult <- runGC(files = c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), settings = TSQXLS.GC, DB = DB, nSlaves = 2)
juist <- to.msp(object = allSamples, file = NULL, settings = metaSetting(settings, "DBconstruction"))

settings2 <- settings


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
  result <- removeDoubleMasses(lapply(pspectra, function(x) cbind(mz = round(allpks[x, 
                                                                                    "mz"], digits = ndigit), allpks[x, c(intensity, "rt")])))
  if (secs2mins) {
    invisible(lapply(result, function(x) cbind(x[, c("mz", 
                                                     intensity)], rt = x[, "rt"]/60)))
  }
  else {
    invisible(result)
  }
}





# Convert MAR xcmsSet to xcms xcmsSet, no method for convertion, since slots are the same, loop through.
xset <- smSet$xcmsSet
for (i in slotNames(xset)) {
  slot(lekkersetjehoor, i) <- slot(xset, i)
}

annxset <- xsAnnotate(xs = lekkersetjehoor, sample = c(1:2), polarity = 'positive')
setjesgroeperen <- groupFWHM(object = annxset)

annxset1 <- xsAnnotate(xs = xcms_peaks, sample = c(1:2), polarity = 'positive')
setjesgroeperen <- groupFWHM(object = annxset)


result3 <- runGC(xset = setjesgroeperen, settings = TSQXLS.GC, DB = MoNa_MSP_metaMS)


function (files, xset, settings, rtrange = NULL, DB = NULL, removeArtefacts = TRUE, 
          findUnknowns = nexp >= mcs, returnXset = FALSE, RIstandards = NULL, 
          nSlaves = 0) 
{
  if (!missing(files)) {
    nexp <- length(files)
  }
  else {
    if (missing(xset)) 
      stop("Either 'files' or 'xset' should be given")
    if (class(xset) == "xcmsSet") 
      stop("xset should be a list of CAMERA-grouped xcmsSet objects, see man page")
    xset.l <- xset
    nexp <- length(xset.l)
  }
  mcs <- max(2, min(metaSetting(settings, "betweenSamples.min.class.size"), 
                    metaSetting(settings, "betweenSamples.min.class.fraction") * 
                      nexp))
  if (findUnknowns & nexp < mcs) {
    stop("Number of samples too small to define unknowns - either provide more samples or change the settings.")
  }
  if (is.null(DB) & !findUnknowns) 
    stop("Nothing to do. Provide a DB or set 'findUnknowns' to TRUE...")
  if (findUnknowns & !is.null(DB) & metaSetting(settings, "betweenSamples.timeComparison") != 
      metaSetting(settings, "match2DB.timeComparison")) 
    stop("Settings error: choose one value for timeComparison only...")
  if (is.null(RIstandards) & ((metaSetting(settings, "betweenSamples.timeComparison") == 
                               "RI" & findUnknowns) | (metaSetting(settings, "match2DB.timeComparison") == 
                                                       "RI" & !is.null(DB)))) 
    stop("Argument RIstandards is mandatory when using RI for matching")
  if (!is.null(RIstandards) & metaSetting(settings, "betweenSamples.timeComparison") == 
      "rt" & metaSetting(settings, "match2DB.timeComparison") == 
      "rt") 
    printWarning("Warning: argument RIstandards provided, but using retention times for matching")
  printString(paste("Experiment of", nexp, "samples"))
  printString(paste("Instrument:", metaSetting(settings, "protocolName")))
  if (length(rtrange) == 2) 
    printString(paste("Retention time range:", rtrange[1], 
                      "to", rtrange[2], "minutes"))
  if (!is.null(DB)) {
    DB.orig <- DB
    DB <- treat.DB(DB.orig)
    printString(paste("Annotation using database of", length(DB), 
                      "spectra"))
  }
  else {
    printString("No annotation performed")
    DB.orig <- NULL
  }
  if (!missing(files)) {
    printString("Performing peak picking and CAMERA")
    xset.l <- peakDetection(files, metaSetting(settings, 
                                               "PeakPicking"), rtrange = rtrange, convert2list = TRUE, 
                            nSlaves = nSlaves)
    allSamples <- lapply(xset.l, runCAMERA, chrom = metaSetting(settings, 
                                                                "chrom"), settings = metaSetting(settings, "CAMERA"))
  }
  else {
    printString("Using xcmsSet object - only doing annotation")
    allSamples <- xset.l
  }
  allSamples.msp <- lapply(allSamples, to.msp, file = NULL, 
                           settings = metaSetting(settings, "DBconstruction"))
  names(allSamples.msp) <- sapply(allSamples, function(x) sampnames(x@xcmsSet))
  if (!is.null(RIstandards)) 
    allSamples.msp <- lapply(allSamples.msp, addRI, RIstandards, 
                             isMSP = FALSE)
  nofeats <- which(sapply(allSamples.msp, length) == 0)
  if ((nnof <- length(nofeats)) > 0) {
    printWarning("Removing", nnof, "injections without any features:\n\t", 
                 paste(names(allSamples)[nofeats], collapse = "\n\t "))
    allSamples.msp <- allSamples.msp[-nofeats]
  }
  allSamples.msp.scaled <- lapply(allSamples.msp, treat.DB, 
                                  isMSP = FALSE)
  if (!is.null(DB)) {
    if (removeArtefacts) {
      printString(paste("Removing artefacts (", paste(metaSetting(settings, 
                                                                  "matchIrrelevants.irrelevantClasses"), collapse = ", "), 
                        ")", sep = ""))
      irrel.idx <- which(sapply(DB, function(x) x$Class) %in% 
                           metaSetting(settings, "matchIrrelevants.irrelevantClasses"))
      if (length(irrel.idx) > 0) {
        subDB <- DB[irrel.idx]
        junkPatterns <- lapply(matchSamples2DB(allSamples.msp.scaled, 
                                               subDB, metaSetting(settings, "matchIrrelevants"), 
                                               quick = TRUE)$annotations, function(x) x[, 
                                                                                        "pattern"])
        allSamples.msp <- mapply(function(x, y) if (length(y) > 
                                                    0) {
          x[-y]
        }
        else {
          x
        }, allSamples.msp, junkPatterns)
        allSamples.msp.scaled <- mapply(function(x, y) if (length(y) > 
                                                           0) {
          x[-y]
        }
        else {
          x
        }, allSamples.msp.scaled, junkPatterns)
        nofeats <- which(sapply(allSamples.msp, length) == 
                           0)
        if ((nnof <- length(nofeats)) > 0) {
          printWarning("Removing", nnof, "injections containing only artefacts:\n\t", 
                       paste(names(allSamples)[nofeats], collapse = "\n\t "))
          allSamples <- allSamples[-nofeats]
        }
        DB <- DB[-irrel.idx]
        DB.orig <- DB.orig[-irrel.idx]
      }
    }
    printString("Matching with database of standards")
    allSam.matches <- matchSamples2DB(allSamples.msp, DB = DB, 
                                      settings = metaSetting(settings, "match2DB"), quick = FALSE)
  }
  else {
    allSam.matches <- NULL
  }
  if (findUnknowns) {
    printString("Matching unknowns across samples")
    allSam.matches <- matchSamples2Samples(allSamples.msp.scaled, 
                                           allSamples.msp, annotations = allSam.matches$annotations, 
                                           settings = metaSetting(settings, "betweenSamples"))
  }
  printString("Formatting results")
  PseudoSpectra <- constructExpPseudoSpectra(allMatches = allSam.matches, 
                                             standardsDB = DB.orig)
  features.df <- getFeatureInfo(stdDB = DB.orig, allMatches = allSam.matches, 
                                sampleList = allSamples.msp)
  ann.df <- getAnnotationMat(exp.msp = allSamples.msp, pspectra = PseudoSpectra, 
                             allMatches = allSam.matches)
  ann.df2 <- sweep(ann.df, 1, sapply(PseudoSpectra, function(x) max(x$pspectrum[, 
                                                                                2])), FUN = "*")
  printString("Done!")
  if (returnXset) {
    list(PeakTable = cbind(data.frame(features.df), data.frame(round(ann.df2))), 
         PseudoSpectra = PseudoSpectra, settings = settings, 
         xset = allSamples, annotation = allSam.matches$annotation, 
         samples.msp = allSamples.msp, SessionInfo = sessionInfo())
  }
  else {
    list(PeakTable = cbind(data.frame(features.df), data.frame(round(ann.df2))), 
         PseudoSpectra = PseudoSpectra, settings = settings, 
         SessionInfo = sessionInfo())
  }
}