# Mummichog for marathon project

data <- readMSData(c('mzxml/Kruid_130/Kruid 130 Zwarte peper 5 191119me_70.mzXML', 'mzxml/Kruid_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML'), mode = 'onDisk')
# data2 <- readMSData('kruiden/Kruid 130 Zwarte peper 5 191119me_70.mzXML', mode = 'onDisk')
param_optimized2 <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0,
                                 snthresh = 100, bw = 2, ppm = 22.35, min_peakwidth = 3, max_peakwidth = 37.5, 
                                 noise = 1000, prefilter = 2, value_of_prefilter = 0.02, minFraction = 0.5, 
                                 minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
                                 family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
                                 mzCenterFun = "wMean")

# param_test <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', mzdiff = 0.01,
#                            snthresh = 10, bw = 2, ppm = 22, min_peakwidth = 5, max_peakwidth = 30, 
#                            noise = 1000, prefilter = 3, value_of_prefilter = 0.02, minFraction = 0.5, 
#                            minSamples = 1, maxFeatures = 100, extra = 1, span = 0.25, smooth = 'loess', 
#                            family = 'gaussian', verbose.columns = FALSE, fitgauss = FALSE, integrate = 1, 
#                            mzCenterFun = "wMean")

smSet <- PerformPeakPicking(data, param = updateRawSpectraParam(param_optimized2))
smSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- smSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
smSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
smSet <- PerformPeakAlignment(smSet, param = updateRawSpectraParam(param_optimized2))
smSet <- PerformPeakFiling(smSet, param = updateRawSpectraParam(param_optimized2))
smSet$xcmsSet@phenoData$sample_group <- gsub(' ', x = smSet$xcmsSet@phenoData$sample_name, replacement = '_')

annParams <- SetAnnotationParam(polarity = 'negative')
annotPeaks <- PerformPeakAnnotation(mSet = smSet, annotaParam = annParams)


maPeaks <- FormatPeakList(annotPeaks = annotPeaks, annParams = annParams, filtIso = F, filtAdducts = F, missPercent = 0.5, includeRT = T)


smSet <- InitDataObjects("pktable", 'stat', paired = FALSE)
smSet <- Read.TextData(mSetObj = smSet, filePath = 'metaboanalyst_input.csv', format = 'colu', lbl.type = 'disc')



kmSet1 <- SanityCheckData(kmSet)

mSetObj <- MetaboAnalystR:::.get.mSet(kmSet1)
mSetObj <- kmSet

function (mSetObj = NA, filePath, format = "rowu", lbl.type = "disc") 
{
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  mSetObj$dataSet$cls.type <- lbl.type
  mSetObj$dataSet$format <- format
  dat <- .readDataTable(filePath)
  if (class(dat) == "try-error" || ncol(dat) == 1) {
    AddErrMsg("Data format error. Failed to read in the data!")
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv) format!")
    AddErrMsg("Please also check the followings: ")
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.")
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.")
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.")
    AddErrMsg("Missing values should be blank or NA without quote.")
    AddErrMsg("Make sure the file delimeters are commas.")
    return(0)
  }
  msg <- NULL
  if (substring(format, 4, 5) == "ts") {
    if (substring(format, 1, 3) == "row") {
      msg <- c(msg, "Samples are in rows and features in columns")
      smpl.nms <- dat[, 1]
      all.nms <- colnames(dat)
      facA.lbl <- all.nms[2]
      cls.lbl <- facA <- dat[, 2]
      facB.lbl <- all.nms[3]
      facB <- dat[, 3]
      conc <- dat[, -c(1:3)]
      var.nms <- colnames(conc)
    }
    else {
      msg <- c(msg, "Samples are in columns and features in rows.")
      all.nms <- dat[, 1]
      facA.lbl <- all.nms[1]
      cls.lbl <- facA <- dat[1, -1]
      facB.lbl <- all.nms[2]
      facB <- dat[2, -1]
      var.nms <- dat[-c(1:2), 1]
      conc <- t(dat[-c(1:2), -1])
      smpl.nms <- rownames(conc)
    }
    if (mSetObj$dataSet$design.type == "time" | mSetObj$dataSet$design.type == 
        "time0") {
      if (!(tolower(facA.lbl) == "time" | tolower(facB.lbl) == 
            "time")) {
        AddErrMsg("No time points found in your data")
        AddErrMsg("The time points group must be labeled as <b>Time</b>")
        return(0)
      }
    }
  }
  else {
    if (substring(format, 1, 3) == "row") {
      msg <- c(msg, "Samples are in rows and features in columns")
      smpl.nms <- dat[, 1]
      dat[, 1] <- NULL
      if (lbl.type == "qc") {
        rownames(dat) <- smpl.nms
        mSetObj$dataSet$orig <- dat
        mSetObj$dataSet$cmpd <- colnames(dat)
        return(1)
      }
      cls.lbl <- dat[, 1]
      conc <- dat[, -1, drop = FALSE]
      var.nms <- colnames(conc)
    }
    else {
      msg <- c(msg, "Samples are in columns and features in rows.")
      var.nms <- dat[-1, 1]
      dat[, 1] <- NULL
      smpl.nms <- colnames(dat)
      cls.lbl <- dat[1, ]
      conc <- t(dat[-1, ])
    }
  }
  mSetObj$dataSet$type.cls.lbl <- class(cls.lbl)
  dat <- NULL
  msg <- c(msg, "The uploaded file is in comma separated values (.csv) format.")
  empty.inx <- is.na(smpl.nms) | smpl.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), 
                        "empty rows</font> were detected and excluded from your data."))
    smpl.nms <- smpl.nms[!empty.inx]
    cls.lbl <- cls.lbl[!empty.inx]
    conc <- conc[!empty.inx, ]
  }
  empty.inx <- is.na(cls.lbl) | cls.lbl == ""
  if (sum(empty.inx) > 0) {
    if (mSetObj$analSet$type != "roc") {
      msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), 
                          "empty labels</font> were detected and excluded from your data."))
      smpl.nms <- smpl.nms[!empty.inx]
      cls.lbl <- cls.lbl[!empty.inx]
      conc <- conc[!empty.inx, ]
    }
    else {
      cls.lbl[is.na(cls.lbl)] <- ""
      msg <- c(msg, paste("<font color=\"orange\">", sum(empty.inx), 
                          "new samples</font> were detected from your data."))
    }
  }
  if (anal.type == "roc") {
    if (length(unique(cls.lbl[!empty.inx])) > 2) {
      AddErrMsg("ROC analysis is only defined for two-group comparisions!")
      return(0)
    }
  }
  if (length(unique(smpl.nms)) != length(smpl.nms)) {
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse = " ")
    AddErrMsg("Duplicate sample names are not allowed!")
    AddErrMsg(dup.nm)
    return(0)
  }
  empty.inx <- is.na(var.nms) | var.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), 
                        "empty features</font> were detected and excluded from your data."))
    var.nms <- var.nms[!empty.inx]
    conc <- conc[, !empty.inx]
  }
  if (length(unique(var.nms)) != length(var.nms)) {
    dup.nm <- paste(var.nms[duplicated(var.nms)], collapse = " ")
    AddErrMsg("Duplicate feature names are not allowed!")
    AddErrMsg(dup.nm)
    return(0)
  }
  if (sum(is.na(iconv(smpl.nms))) > 0) {
    na.inx <- is.na(iconv(smpl.nms))
    nms <- paste(smpl.nms[na.inx], collapse = "; ")
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", 
                    nms, collapse = " "))
    return(0)
  }
  if (sum(is.na(iconv(var.nms))) > 0) {
    na.inx <- is.na(iconv(var.nms))
    nms <- paste(var.nms[na.inx], collapse = "; ")
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", 
                    nms, collapse = " "))
    return(0)
  }
  smpl.nms <- MetaboAnalystR:::CleanNames(smpl.nms, "sample_name")
  orig.var.nms <- var.nms
  var.nms <- MetaboAnalystR:::CleanNames(var.nms, "var_name")
  names(orig.var.nms) <- var.nms
  cls.lbl <- MetaboAnalystR:::ClearStrings(as.vector(cls.lbl))
  rownames(conc) <- smpl.nms
  colnames(conc) <- var.nms
  if (mSetObj$dataSet$paired) {
    mSetObj$dataSet$orig.cls <- mSetObj$dataSet$pairs <- cls.lbl
  }
  else {
    if (lbl.type == "disc") {
      if (min(table(cls.lbl)) < 3) {
        AddErrMsg(paste("A total of", length(levels(as.factor(cls.lbl))), 
                        "groups found with", length(smpl.nms), "samples."))
        AddErrMsg("At least three replicates are required in each group!")
        AddErrMsg("Or maybe you forgot to specify the data format?")
        return(0)
      }
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.factor(as.character(cls.lbl))
      if (substring(format, 4, 5) == "ts") {
        mSetObj$dataSet$facA.type <- is.numeric(facA)
        mSetObj$dataSet$orig.facA <- mSetObj$dataSet$facA <- as.factor(as.character(facA))
        mSetObj$dataSet$facA.lbl <- facA.lbl
        mSetObj$dataSet$facB.type <- is.numeric(facB)
        mSetObj$dataSet$orig.facB <- mSetObj$dataSet$facB <- as.factor(as.character(facB))
        mSetObj$dataSet$facB.lbl <- facB.lbl
      }
    }
    else {
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- tryCatch({
        as.numeric(cls.lbl)
      }, warning = function(na) {
        print("Class labels must be numeric and continuous!")
        return(0)
      })
      if (mSetObj$dataSet$cls == 0) {
        AddErrMsg("Class labels must be numeric and continuous!")
        return(0)
      }
    }
  }
  if (mSetObj$dataSet$type == "conc") {
    mSetObj$dataSet$cmpd <- var.nms
  }
  mSetObj$dataSet$mumType <- "table"
  mSetObj$dataSet$orig.var.nms <- orig.var.nms
  mSetObj$dataSet$orig <- conc
  mSetObj$msgSet$read.msg <- c(msg, paste("The uploaded data file contains ", 
                                          nrow(conc), " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSetObj$dataSet$type)), 
                                          ") data matrix.", sep = ""))
  return(.set.mSet(mSetObj))
}