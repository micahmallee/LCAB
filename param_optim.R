raw_data1 <- ImportRawMSData(foldername = 'C://Users/Micah/Documents/Data_IBD/', mode = 'inMemory', plotSettings = SetPlotParam(Plot = F))
raw_data2 <- readMSData(files = 'mzxml/KRUID_131/Kruid 131 Zwarte peper 6 191119me_71.mzXML', mode = 'inMemory', msLevel. = 1)
raw_data = raw_data2

Sys.setenv('R_REMOTES_NO_ERRORS_FROM_WARNINGS' = 'true')

POtest <- PerformParamsOptimization(raw_data = raw_data, param = SetPeakParam(platform = 'general', Peak_method = 'centWave'), ncore = 1)

raw_data <- PerformDataTrimming('c://Users/Micah/Documents/QC_IBD/',rt.idx = 0.2)
param_initial <- SetPeakParam(platform = "UPLC-Q/E") 
param_optimized <- PerformParamsOptimization(raw_data, param = param_initial, ncore = 1)

param_optimized <-OptiLCMS::PerformParamsOptimization(raw_data, param = OptiLCMS::SetPeakParam(), ncore = 1)


function (raw_data, param = p0, method = "DoE", ncore = 4)
{
  start.time <- Sys.time()
  if (missing(param)) {
    stop("Please provide the param with 'SetPeakParam' function !")
  }
  else if (missing(raw_data)) {
    stop("Please provide the data of MSnExp format!")
  }
  else if (missing(method)) {
    method <- "DoE"
    print("DoE Optimization Starting Now!")
  }
  if (missing(ncore)) {
    ncore <- detectCores()
    message("'ncore' is absent, will use 2/3 CPU threads of total!")
  }
  if (ncore == 1) {
    if (.Platform$OS.type == "unix") {
      register(bpstart(MulticoreParam(ncore)))
    }
    else {
      register(bpstart(SnowParam(ncore)))
    }
  }
  if (param[["Peak_method"]] == "centWave") {
    print("Evaluating Noise level...")
    p2 <- Noise_evaluate(raw_data)
    param[["ppm"]] <- round(p2$ppm, 2)
    param[["noise"]] <- round(p2$noise, 2)
    param[["prefilter"]] <- round(p2$prefilter, 2)
    param[["value_of_prefilter"]] <- round(p2$value_of_prefilter, 
                                           2)
  }
  if (method == "DoE") {
    p1 <- MetaboAnalystR:::optimize.xcms.doe(raw_data, param = param, ncore = ncore)
  }
  if (method == "OVAT") {
    stop("Only DoE is supported for now. Other Optimization Model will be supported later.")
    p1 <- optimize.xcms.ovat(raw_data, param = param, ncore = ncore)
  }
  end.time <- Sys.time()
  message("Time Spent In Total:", round((as.numeric(end.time) - 
                                           as.numeric(start.time))/60, 1), "mins", "\n")
  return(p1)
}








function (raw_data, param, ncore = 8) 
{
  if (is.null(param)) {
    stop("Please provide the param with 'SetPeakParam' function !")
  }
  else {
    Parameters <- param
  }
  if (Parameters$Peak_method == "centWave" && Parameters$RT_method == 
      "peakgroup") {
    Parameters$value_of_prefilter <- Parameters$value_of_prefilter
    Parameters$max_peakwidth <- c(Parameters$max_peakwidth * 
                                    0.5, Parameters$max_peakwidth * 2)
    Parameters$min_peakwidth <- c((Parameters$min_peakwidth) * 
                                    0.5, (Parameters$min_peakwidth) * 2)
    Parameters$mzdiff <- c(-Parameters$mzdiff * 1.2, Parameters$mzdiff * 
                             1.2)
    Parameters$snthresh <- c(Parameters$snthresh * 0.75, 
                             Parameters$snthresh * 5)
    Parameters$bw <- c(Parameters$bw * 0.5, Parameters$bw * 
                         1.5)
  }
  result <- MetaboAnalystR:::optimizxcms.doe.peakpicking(object = raw_data, 
                                        params = Parameters, BPPARAM = bpparam(), nSlaves = ncore, 
                                        subdir = NULL, plot = F)
  optimizedxcmsObject <- result$best_settings$xset
  message("Optimization Finished !")
  peakParams2 <- list()
  peakParams2$best_parameters <- result[["best_settings"]][["parameters"]]
  peakParams2$data <- optimizedxcmsObject
  message("Parameters Optimization Finished !")
  return(peakParams2)
}










function (object = NULL, params = params, BPPARAM = bpparam(), 
          nSlaves = 4, plot = F, ...) 
{
  isotopeIdentification = c("IPO")
  centWave <- is.null(params$fwhm)
  history <- list()
  iterator = 1
  best_range <- 0.25
  object_mslevel <- MetaboAnalystR:::PeakPicking_prep(object)
  while (iterator < 20) {
    print(paste0("Round:", iterator))
    message("DoE Running Begin...")
    mSet_OPT <- ExperimentsCluster_doe(object = object, object_mslevel = object_mslevel, 
                                       params = params, isotopeIdentification = isotopeIdentification, 
                                       BPPARAM = BPPARAM, nSlaves = nSlaves)
    PPS.set <- as.numeric(sapply(1:nrow(mSet_OPT[["response"]]), 
                                 FUN = function(x) {
                                   mSet_OPT[["response"]][x, 5]
                                 }))
    CV.set <- as.numeric(sapply(1:nrow(mSet_OPT[["response"]]), 
                                FUN = function(x) {
                                  mSet_OPT[["response"]][x, 6]
                                }))
    RCS.set <- as.numeric(sapply(1:nrow(mSet_OPT[["response"]]), 
                                 FUN = function(x) {
                                   mSet_OPT[["response"]][x, 7]
                                 }))
    GS.set <- as.numeric(sapply(1:nrow(mSet_OPT[["response"]]), 
                                FUN = function(x) {
                                  mSet_OPT[["response"]][x, 8]
                                }))
    GaussianSI.set <- as.numeric(sapply(1:nrow(mSet_OPT[["response"]]), 
                                        FUN = function(x) {
                                          mSet_OPT[["response"]][x, 9]
                                        }))
    index.set <- list(CV = CV.set, RCS = RCS.set, GS = GS.set, 
                      GaussianSI = GaussianSI.set)
    CV.set.normalized <- (CV.set - min(CV.set))/(max(CV.set) - 
                                                   min(CV.set))
    RCS.set.normalized <- (RCS.set - min(RCS.set))/(max(RCS.set) - 
                                                      min(RCS.set))
    GS.set.normalized <- (GS.set - min(GS.set))/(max(GS.set) - 
                                                   min(GS.set))
    GaussianSI.set.normalized <- GaussianSI.set
    QCoE <- CV.set.normalized * 0.2 + RCS.set.normalized * 
      0.4 + GS.set.normalized * 0.4
    QS <- PPS.set * QCoE * GaussianSI.set.normalized
    tmp_matrix <- mSet_OPT[["response"]]
    tmp_matrix <- cbind(tmp_matrix, PPS.set, CV.set.normalized, 
                        RCS.set.normalized, GS.set.normalized, GaussianSI.set.normalized, 
                        QCoE, QS)
    mSet_OPT[["response"]] <- tmp_matrix
    message("Round ", iterator, " Finished !")
    mSet_OPT <- Statistic_doe(object = object, object_mslevel = object_mslevel, 
                              isotopeIdentification = isotopeIdentification, BPPARAM = BPPARAM, 
                              subdir = NULL, plot = F, mSet_OPT = mSet_OPT, iterator = iterator, 
                              index.set = index.set, useNoise = params[["noise"]])
    history[[iterator]] <- mSet_OPT
    params <- mSet_OPT$params
    if (!resultIncreased_doe(history)) {
      message("No Increase Stopping !")
      maxima <- 0
      max_index <- 1
      for (i in 1:length(history)) {
        if (history[[i]]$max_settings[1] > maxima) {
          maxima <- history[[i]]$max_settings[1]
          max_index <- i
        }
      }
      xcms_parameters <- as.list(decodeAll(history[[max_index]]$max_settings[-1], 
                                           history[[max_index]]$params$to_optimize))
      xcms_parameters <- combineParams(xcms_parameters, 
                                       params$no_optimization)
      if (!is.list(xcms_parameters)) 
        xcms_parameters <- as.list(xcms_parameters)
      best_settings <- list()
      best_settings$parameters <- xcms_parameters
      best_settings$xset <- history[[max_index]]$xset
      target_value <- history[[max_index]]$QS
      best_settings$result <- target_value
      history$best_settings <- best_settings
      message("best parameter settings:")
      message(paste(rbind(paste(names(xcms_parameters), 
                                sep = "", ": "), paste(xcms_parameters, 
                                                       sep = "", "\n")), sep = ""))
      return(history)
    }
    for (i in 1:length(params$to_optimize)) {
      parameter_setting <- mSet_OPT$max_settings[i + 1]
      bounds <- params$to_optimize[[i]]
      fact <- names(params$to_optimize)[i]
      min_factor <- ifelse(fact == "min_peakwidth", 
                           3, ifelse(fact == "mzdiff", ifelse(centWave, 
                                                              -0.02, 0.001), ifelse(fact == "step", 
                                                                                    5e-04, ifelse(fact == "bw", 2, ifelse(fact == 
                                                                                                                            "snthresh", 5, 1)))))
      step_factor <- ifelse(is.na(parameter_setting), 1.2, 
                            ifelse((abs(parameter_setting) < best_range), 
                                   0.8, ifelse(parameter_setting == -1 & decode(-1, 
                                                                                params$to_optimize[[i]]) == min_factor, 0.8, 
                                               1)))
      step <- (diff(bounds)/2) * step_factor
      if (is.na(parameter_setting)) 
        parameter_setting <- 0
      new_center <- decode(parameter_setting, bounds)
      if ((new_center - min_factor) > step) {
        new_bounds <- c(new_center - step, new_center + 
                          step)
      }
      else {
        new_bounds <- c(min_factor, 2 * step + min_factor)
      }
      names(new_bounds) <- NULL
      if (names(params$to_optimize)[i] == "steps" | 
          names(params$to_optimize)[i] == "prefilter") {
        params$to_optimize[[i]] <- round(new_bounds, 
                                         0)
      }
      else {
        params$to_optimize[[i]] <- new_bounds
      }
    }
    if (centWave) {
      if (!is.null(params$to_optimize$min_peakwidth) | 
          !is.null(params$to_optimize$max_peakwidth)) {
        pw_min <- ifelse(is.null(params$to_optimize$min_peakwidth), 
                         params$no_optimization$min_peakwidth, max(params$to_optimize$min_peakwidth))
        pw_max <- ifelse(is.null(params$to_optimize$max_peakwidth), 
                         params$no_optimization$max_peakwidth, min(params$to_optimize$max_peakwidth))
        if (pw_min >= pw_max) {
          additional <- abs(pw_min - pw_max) + 1
          if (!is.null(params$to_optimize$max_peakwidth)) {
            params$to_optimize$max_peakwidth <- params$to_optimize$max_peakwidth + 
              additional
          }
          else {
            params$no_optimization$max_peakwidth <- params$no_optimization$max_peakwidth + 
              additional
          }
        }
      }
    }
    params <- attachList(params$to_optimize, params$no_optimization)
    iterator <- iterator + 1
  }
  params <- attachList(params$to_optimize, params$no_optimization)
  return(history)
}


{
  centroided <- all(centroided(object))
  if (is.na(centroided)) {
    suppressWarnings(centroided <- isCentroided(object[[ceiling(length(object)/3)]]))
  }
  if (is.na(centroided) || !centroided) 
    warning("Your data appears to be not centroided! CentWave", 
            " works best on data in centroid mode.")
  object_mslevel_i <- splitByFile(object = object, f = factor(c(1:length(object@phenoData@data[["sample_name"]]))))
  object_mslevel_o <- bplapply(object_mslevel_i, FUN = function(x) {
    x@assayData
  })
  message("Data Spliting Finished !")
  message("Peak Preparing Begin...")
  object_mslevel_name <- bplapply(object_mslevel_o, ls)
  object_mslevell <- object_mslevel <- list()
  for (j in 1:length(object_mslevel_o)) {
    for (i in 1:length(object_mslevel_name[[j]])) {
      object_mslevell[[i]] <- object_mslevel_o[[j]][[object_mslevel_name[[j]][i]]]
    }
    names(object_mslevell) <- object_mslevel_name[[j]]
    object_mslevel[[j]] <- object_mslevell
  }
  for (i in 1:length(object_mslevel)) {
    ncount <- as.numeric(which(sapply(names(object_mslevel[[i]]), 
                                      is.na)))
    if (!identical(ncount, numeric(0))) {
      object_mslevel[[i]] <- object_mslevel[[i]][-ncount]
    }
  }
  message("Peak Preparing Done !")
  return(object_mslevel)
}








