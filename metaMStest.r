# Testing MetaMS
library(metaMS)
library(metaMSdata)
data(threeStdsDB)
data(FEMsettings)

cdfdir <- system.file('extdata', package = 'metaMSdata')
cdffiles <- list.files(cdfdir, pattern = "_GC_", full.names = T, ignore.case = T)
result <- runGC(files = cdffiles, settings = TSQXLS.GC, DB = DB, nSlaves = 2)

<<<<<<< HEAD
testresult <- runGC(files = c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), settings = TSQXLS.GC, DB = MoNA_DB, nSlaves = 2)
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


# Peak picking
GCset <- peakDetection(files = cdffiles, settings = metaSetting(TSQXLS.GC, 'PeakPicking'), convert2list = T, nSlaves = 2)
# Definition of pseudospectra
allSamples <- lapply(GCset, runCAMERA, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA"))
allSamples.msp <- lapply(allSamples, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))
sapply(allSamples.msp, length)
plotPseudoSpectrum(allSamples.msp[[1]][[26]])
# Annotation
DB.treated <- treat.DB(DB)
allSam.matches <- matchSamples2DB(allSamples.msp, DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)
allSam.matches
matchExpSpec(allSamples.msp[[1]][[4]], DB.treated, DB.treated = TRUE, plotIt = TRUE)
# Unknowns
allSamples.msp.scaled <- lapply(allSamples.msp, treat.DB, isMSP = F)
allSamples.matches <- matchSamples2Samples(allSamples.msp.scaled, allSamples.msp, annotations = allSam.matches$annotations, settings = metaSetting(TSQXLS.GC, "betweenSamples"))
names(allSamples.matches)
# Output
features.df <- getFeatureInfo(stdDB = DB, allMatches = allSamples.matches, sampleList = allSamples.msp)
features.df[, c(1:3, ncol(features.df) - 2:0)]
PseudoSpectra <- constructExpPseudoSpectra(allMatches = allSamples.matches, standardsDB = DB)
ann.df <- getAnnotationMat(exp.msp = allSamples.msp, pspectra = PseudoSpectra, allMatches = allSamples.matches)
ann.df
ann.df2 <- sweep(ann.df, 1, sapply(PseudoSpectra, function(x) max(x$pspectrum[, 2])), FUN = "*")
ann.df2


oke <- xsAnnotate(xs = xset)
oke1 <- groupFWHM(object = oke)

xsetCam <- runCAMERA(xset, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA"))

xsetmsp <- to.msp(xsetgrouped, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))


# Convert MAR xcmsSet to xcms xcmsSet, no method for conversion, since slots are the same, loop through.
xset <- smSet$xcmsSet
annxset <- xsAnnotate(xs = xset, sample = c(1:2), polarity = 'positive')
setjesgroeperen <- groupFWHM(object = annxset)

# for (i in slotNames(xset)) {
#   slot(lekkersetjehoor, i) <- slot(xset, i)
# }

result3 <- runGC(xset = setjesgroeperen, settings = TSQXLS.GC, DB = MoNA_DB)




# GroupFWHM
function (object, sigma = 6, perfwhm = 0.6, intval = "maxo") 
{
  if (!class(object) == "xsAnnotate") {
    stop("no xsAnnotate object")
  }
  if (!sum(intval == c("into", "intb", "maxo"))) {
    stop("unknown intensity value!")
  }
  sample <- object@sample
  pspectra <- list()
  psSamples <- NA
  cat("Start grouping after retention time.\n")
  if (object@groupInfo[1, "rt"] == -1) {
    warning("Warning: no retention times avaiable. Do nothing\n")
    return(invisible(object))
  }
  else {
    if (is.na(sample[1]) || length(object@xcmsSet@filepaths) > 
        1) {
      if (is.na(sample[1])) {
        index <- 1:length(object@xcmsSet@filepaths)
      }
      else {
        index <- sample
      }
      gvals <- groupval(object@xcmsSet)[, index, drop = FALSE]
      peakmat <- object@xcmsSet@peaks
      groupmat <- groups(object@xcmsSet)
      maxo <- as.numeric(apply(gvals, 1, function(x, peakmat) {
        val <- na.omit(peakmat[x, intval])
        if (length(val) == 0) {
          return(NA)
        } else {
          return(max(val))
        }
      }, peakmat))
      maxo[which(is.na(maxo))] <- -1
      maxo <- cbind(1:length(maxo), maxo)
      int.max <- as.numeric(apply(gvals, 1, function(x, 
                                                     peakmat) {
        which.max(peakmat[x, intval])
      }, peakmat))
      peakrange <- matrix(apply(gvals, 1, function(x, peakmat) {
        val <- peakmat[x, intval]
        if (length(na.omit(val)) == 0) {
          return(c(0, 1))
        }
        else {
          return(peakmat[x[which.max(val)], c("rtmin", 
                                              "rtmax")])
        }
      }, peakmat), ncol = 2, byrow = TRUE)
      colnames(peakrange) <- c("rtmin", "rtmax")
      while (length(maxo) > 0) {
        iint <- which.max(maxo[, 2])
        rtmed <- groupmat[iint, "rtmed"]
        rt.min <- peakrange[iint, "rtmin"]
        rt.max <- peakrange[iint, "rtmax"]
        hwhm <- ((rt.max - rt.min)/sigma * 2.35 * perfwhm)/2
        irt <- which(groupmat[, "rtmed"] > (rtmed - hwhm) & 
                       groupmat[, "rtmed"] < (rtmed + hwhm))
        if (length(irt) > 0) {
          idx <- maxo[irt, 1]
          pspectra[[length(pspectra) + 1]] <- idx
          psSamples[length(pspectra)] <- index[int.max[maxo[iint, 
                                                            1]]]
          maxo <- maxo[-irt, , drop = FALSE]
          groupmat <- groupmat[-irt, , drop = FALSE]
          peakrange <- peakrange[-irt, , drop = FALSE]
        }
        else {
          idx <- maxo[iint, 1]
          cat("Warning: Feature ", idx, " looks odd for at least one peak. Please check afterwards.\n")
          pspectra[[length(pspectra) + 1]] <- idx
          psSamples[length(pspectra)] <- index[int.max[maxo[iint, 
                                                            1]]]
          maxo <- maxo[-iint, , drop = FALSE]
          groupmat <- groupmat[-iint, , drop = FALSE]
          peakrange <- peakrange[-iint, , drop = FALSE]
        }
      }
    }
    else {
      peakmat <- object@xcmsSet@peaks
      maxo <- peakmat[, intval]
      maxo <- cbind(1:length(maxo), maxo)
      while (length(maxo) > 0) {
        iint <- which.max(maxo[, 2])
        rtmed <- peakmat[iint, "rt"]
        rt.min <- peakmat[iint, "rtmin"]
        rt.max <- peakmat[iint, "rtmax"]
        hwhm <- ((rt.max - rt.min)/sigma * 2.35 * perfwhm)/2
        irt <- which(peakmat[, "rt"] > (rtmed - hwhm) & 
                       peakmat[, "rt"] < (rtmed + hwhm))
        if (length(irt) > 0) {
          idx <- maxo[irt, 1]
          pspectra[[length(pspectra) + 1]] <- idx
          maxo <- maxo[-irt, , drop = FALSE]
          peakmat <- peakmat[-irt, , drop = FALSE]
        }
        else {
          idx <- maxo[iint, 1]
          cat("Warning: Feature ", idx, " looks odd for at least one peak. Please check afterwards.\n")
          pspectra[[length(pspectra) + 1]] <- idx
          maxo <- maxo[-iint, , drop = FALSE]
          peakmat <- peakmat[-iint, , drop = FALSE]
        }
      }
      psSamples <- rep(sample, length(pspectra))
    }
    object@pspectra <- pspectra
    object@psSamples <- psSamples
    cat("Created", length(object@pspectra), "pseudospectra.\n")
  }
  return(invisible(object))
}


=======
testresult <- runGC(files = c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), settings = TSQXLS.GC, DB = DB, nSlaves = 2)
juist <- to.msp(object = allSamples, file = NULL, settings = metaSetting(settings, "DBconstruction"))

settings2 <- settings
>>>>>>> 57efe2f095d2926b7a494efe3b7d235598ea8ead


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