# Testing MetaMS
library(metaMS)
library(metaMSdata)
data(threeStdsDB)
data(FEMsettings)

# Complete metaMS workflow
cdfdir <- system.file('extdata', package = 'metaMSdata')
cdffiles <- list.files(cdfdir, pattern = "_GC_", full.names = T, ignore.case = T)
# result1 <- runGC(files = cdffiles, settings = TSQXLS.GC, DB = mmsdb, nSlaves = 2)

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
DB.treated <- treat.DB(mmsdb)
allSam.matches <- matchSamples2DB(allSamples.msp, DB = DB.treated2, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)
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



sampletje <- runCAMERA(xset = GCset[[1]], chrom = 'GC', settings = metaSetting(TSQXLS.GC, "CAMERA"))
mspsampletje <- to.msp(object = sampletje, settings = metaSetting(TSQXLS.GC, 'DBconstruction'))
DB.treated <- treat.DB(DB)
sampletjematches <- matchSamples2DB(list(mspsampletje), DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)




# Own data:
mmsdb <- read.msp(file = 'MoNA-export-GC-MS_Spectra.msp')
mmssettings <- TSQXLS.GC
## Make data suitable for runCAMERA? or for to.msp # xsannotat + groupFWHM = runCAMERA
xsetxcms <- xsAnnotate(xs = smSet$xcmsSet)
### Create pseudospectra
groupedxset <- groupFWHM(object = xsetxcms)
### Correct format peakinfo
convertedgroupedxset <- groupedxset
convertedgroupedxset@groupInfo <- convertedgroupedxset@xcmsSet@peaks
### Convert pseudospectra to MSP format
xsetmsp <- to.msp(object = convertedgroupedxset, settings = metaSetting(TSQXLS.GC, 'DBconstruction'))
length(xsetmsp)
### Annotation
DB.treated <- treat.DB(DB)
nestedxsetmsp <- list(xsetmsp)
xsetallmatches <- matchSamples2DB(xset.msp = nestedxsetmsp, DB = DB.treated, settings = metaSetting(TSQXLS.GC, "match2DB"), quick = FALSE)
xsetallmatches <- matchSamples2DB(xset.msp = nestedxsetmsp, DB = DB.treated, settings = settings, quick = FALSE)
xsetallmatches
matchExpSpec(xsetallmatches[[1]][[1]], DB.treated, DB.treated = TRUE, plotIt = TRUE)
## Unknowns
nestedxsetmsp.scaled <- lapply(nestedxsetmsp, treat.DB, isMSP = F)
nestedxsetmsp.scaled <- list(nestedxsetmsp.scaled)

####################################################################
xsetallmatches_unknowns <- matchSamples2Samples(nestedxsetmsp.scaled, nestedxsetmsp, annotations = xsetallmatches$annotations, settings = metaSetting(TSQXLS.GC, "betweenSamples"))
names(xsetallmatches_unknowns)


# Convert MAR xcmsSet to xcms xcmsSet, no method for conversion, since slots are the same, loop through.
for (i in slotNames(xset)) {
  slot(lekkersetjehoor, i) <- slot(xset, i)
}









xset.msp <- allSamples.msp
DB = DB.treated2
settings = metaSetting(TSQXLS.GC, 'match2DB')
quick = F


matchSamples2DB <- function(xset.msp, DB, settings, quick) {
  if (settings$timeComparison == "RI") {
    standard.rts <- sapply(DB, function(x) x$std.RI)
    ## rt.matches is a list of two-column matrices, one for each xset.msp element. The first column gives the DB entry that
    ## matches with the sample entry in the second column. This is very fast to calculate.
    rt.matches <- lapply(1:length(xset.msp), function(ii) {
      group.rts <- sapply(xset.msp[[ii]], function(x) mean(x[, "RI"]))
      which(abs(outer(standard.rts, group.rts, FUN = "-")) < settings$RIdiff, arr.ind = TRUE)
    })
  } else {
    standard.rts <- sapply(DB, function(x) x$std.rt)
    
    rt.matches <- lapply(1:length(xset.msp), function(ii) {
      group.rts <- sapply(xset.msp[[ii]], function(x) mean(x[, "rt"]))
      which(abs(outer(standard.rts, group.rts, FUN = "-")) < settings$rtdiff, arr.ind = TRUE)
    })
  }
  
  if (quick) {
    ## xset.msp has already been scaled
    match.results <- lapply(1:length(xset.msp), function(ii) {
      if (nrow(rt.matches[[ii]] > 0)) {
        result <- matrix(0, length(DB), length(xset.msp[[ii]]))
        for (i in 1:nrow(rt.matches[[ii]])) {
          DB.idx <- rt.matches[[ii]][i, 1]
          sample.idx <- rt.matches[[ii]][i, 2]
          result[DB.idx, sample.idx] <- mzmatch(DB[[DB.idx]]$pspectrum, xset.msp[[ii]][[sample.idx]])
        }
      }
      result
    })
  } else {
    ## scaling is done for each comparison separately, since high mz values may be removed depending on MonoMW of the standard
    ## compound. Slower, obviously.
    match.results <- lapply(1:length(xset.msp), function(ii) {
      result <- matrix(0, length(DB), length(xset.msp[[ii]]))
      if (nrow(rt.matches[[ii]] > 0)) {
        for (i in 1:nrow(rt.matches[[ii]])) {
          DB.idx <- rt.matches[[ii]][i, 1]
          sample.idx <- rt.matches[[ii]][i, 2]
          exp.pat <- xset.msp[[ii]][[sample.idx]]
          
          MWlimit <- DB[[DB.idx]]$monoMW + 4
          ## sometimes, with manually added spectra, monoMW is not present... then we use everything up until the highest mass present.
          if (length(MWlimit) == 0) 
            MWlimit <- max(DB[[DB.idx]]$pspectrum[, 1])
          
          ok.mz <- which(exp.pat[, "mz"] <= MWlimit)
          if (length(ok.mz) > settings$minfeat) {
            exp.pat <- treat.DB(list(exp.pat[ok.mz, ]), isMSP = FALSE)
            result[DB.idx, sample.idx] <- metaMS:::mzmatch(DB[[DB.idx]]$pspectrum, exp.pat[[1]])
          }
        }
      }
      result
    })
  }
  names(match.results) <- names(xset.msp)
  
  annotations <- lapply(match.results, function(xx) {
    sapply(1:ncol(xx), function(ii) which(xx[, ii] > settings$simthresh))
  })
  
  list(annotations = mapply(annotations2tab, annotations, match.results, SIMPLIFY = FALSE))
}



allsimilarities <-  lapply(1:length(xsetmsp), function(x) {
  for (i in 1:length(DB.treated2)) {
    SpectrumSimilarity(spec.top = xsetmsp[[x]][, 1:2], spec.bottom = DB.treated2[[i]]$pspectrum, print.graphic = F, print.alignment = F)
  }
})



