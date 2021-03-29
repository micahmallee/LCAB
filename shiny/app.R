

library(shiny)
library(MetaboAnalystR)
library(shinydashboard)
library(shinydashboardPlus)
library(xcms)
library(magrittr)
library(shinycssloaders)
library(shinyjs)
library(crosstalk)
library(plotly)
library(metaMS)
library(OrgMassSpecR)
library(splashR)
library(stringr)

# Increase max upload size to 100 MB
options(shiny.maxRequestSize=100*4096^2)

# Load MoNA_DB and SPLASH values
mona_msp <- readRDS(file = 'data/mona_msp')
mona_secondblocks <- readRDS(file = 'data/mona_secondblocks')


# Define UI.
# The front end of the webapp is created here.
ui <- dashboardPagePlus(
  skin = 'midnight',
  header = dashboardHeaderPlus(
    title = tagList(
      tags$span(
        class = "logo-lg", "MetabOracle"
      )
    )
  ),
  # sidenavbar
  sidebar = dashboardSidebar(
    collapsed = T, 
    width = 340,
    sidebarMenu(collapsed = T,
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('dashboard')),
      menuItem('Peak detection', tabName = 'peakpicking', icon = icon('search')),
      menuItem('Clustering', tabName = 'clustering', icon = icon('align-center')),
      menuItem('Metabolite Set Enrichment Analysis', tabName = 'msea', icon = icon('align-center')),
      menuItem('Statistical analysis', tabName = 'statistics', icon = icon('align-center'))
    )
  ),
  # Main content
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem('dashboard', 
              fluidRow(
                column(width = 3,
                  box(width = 12,
                    title = 'Load data',
                    helpText("Welcome to MetabOracle! Please upload your mzXML/mzML/netCDF data below."),
                    fileInput(inputId = 'data_input', label = 'Select your file(s):', multiple = T, accept = c('.mzXML', 'mzxml', '.CDF', '.cdf', '.mzml', 'mzML')),
                    # sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
                    radioButtons(inputId = 'inspect_trim', label = 'Inspect or upload data.',
                                 choices = list('Inspect' = 1,
                                                'Upload' = 2),
                                 selected = 2, inline = T),
                    actionButton(inputId = 'run', label = 'Run')
                  )),
                column(width = 9,
                       box(
                         width = 12,
                         imageOutput(outputId = 'inspect_plot') # %>% withSpinner()
                       ),
                       box(
                         width = 12,
                         title = 'File information:',
                         tableOutput(outputId = 'file')
                       )))),
      tabItem('peakpicking',
              fluidRow(
                column(
                  width = 4,
                  box(id = 'testje', width = 6,
                      title = 'Input parameters for peak picking',
                      numericInput(inputId = 'ppm', label = 'ppm', value = NA, min = 0),
                      bsTooltip(id = 'ppm', title = 'Maximum tolerated fluctuation of m/z value (ppm) from scan to scan - depends on the mass spectrometer accuracy', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'noise', label = 'Noise', value = 1000, min = 0),
                      bsTooltip(id = 'noise', title = 'Each centroid must be greater than the "noise" intensity value', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = 5, min = 0),
                      bsTooltip(id = 'min_peakwidth', title = 'Minimum peak width in seconds', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'mz_diff', label = 'mz diff', value = 0.05, min = -0.01),
                      bsTooltip(id = 'mz_diff', title = 'Minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = 20, min = 0),
                      bsTooltip(id = 'max_peakwidth', title = 'Maximum peak width in seconds', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'snthresh', label = 'Signal to noise threshold', value = 10, min = 0),
                      bsTooltip(id = 'snthresh', title = 'Signal to noise ratio cut-off (intensity)', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'prefilter', label = 'Prefilter', value = 3, min = 0),
                      bsTooltip(id = 'prefilter', title = 'A peak must be present in x scans with an intensity greater than the value of the prefilter', placement = 'left', trigger = 'hover'),
                      numericInput(inputId = 'v_prefilter', label = 'Value of prefilter', value = 100, min = 0),
                      bsTooltip(id = 'v_prefilter', title = 'Value of prefilter', placement = 'left', trigger = 'hover'),
                      actionButton(inputId = 'peakdetectrun', label = 'Perform peak detection'),
                      actionButton(inputId = 'paramdetectrun', label = 'Detect parameters automatically')
                  ),
                  div(id = 'align_param_box',
                      box(width = 6,
                          title = 'Input parameters for peak alignment',
                          selectInput(inputId = 'rtmethod', label = 'Method', choices = list('loess' = 'loess', 'obiwarp' = 'obiwarp'),selected = 'loess'),
                          # loess
                          div(id = 'loessparams',
                              numericInput(inputId = 'extra', label = 'Extra', value = 1, min = 0),
                              bsTooltip(id = 'extra', title = 'Number of "extra" peaks used to define reference peaks (or well-behaved peaks) for modeling time deviation. Number of Peaks > number of samples.', placement = 'left', trigger = 'hover'),
                              numericInput(inputId = 'span', label = 'Span', value = 0.25, min = 0),
                              bsTooltip(id = 'span', title = 'Degree of smoothing of the loess model. 0.2 to 1', placement = 'left', trigger = 'hover')),
                          # obiwarp
                          div(id = 'obiwarpparams',
                              numericInput(inputId = 'prof_step', label = 'profStep', value = 100, min = 0),
                              bsTooltip(id = 'prof_step', title = 'Prof step', placement = 'left', trigger = 'hover')),
                          # Grouping:
                          numericInput(inputId = 'bw', label = 'Bandwith', value = 10, min = 1),
                          bsTooltip(id = 'bw', title = 'Standart deviation of the gaussian metapeak that groups peaks together.', placement = 'left', trigger = 'hover'),
                          numericInput(inputId = 'min_fraction', label = 'minFraction', value = 0.5, min = 0),
                          bsTooltip(id = 'min_fraction', title = 'Minimal fraction of samples a feature has to be present in', placement = 'left', trigger = 'hover'),
                          numericInput(inputId = 'min_samples', label = 'minSamples', value = 1, min = 1),
                          bsTooltip(id = 'min_samples', title = 'Minimum amount of samples a feature has to be present in', placement = 'left', trigger = 'hover'),
                          selectInput(inputId = 'fitgauss', label = 'Fitgauss', choices = list('False' = 'FALSE', 'True' = 'TRUE'),selected = 'FALSE'),
                          bsTooltip(id = 'fitgauss', title = 'If true, a gaussian is fitted to each peak', placement = 'left', trigger = 'hover'),
                          selectInput(inputId = 'verbose_columns', label = 'Verbose columns', choices = list('False' = 'FALSE', 'True' = 'TRUE'),selected = 'FALSE'),
                          bsTooltip(id = 'verbose_columns', title = 'If true additional peak meta data columns are returned', placement = 'left', trigger = 'hover'),
                          numericInput(inputId = 'integrate', label = 'Integrate', value = 1, min = 0, max = 1),
                          bsTooltip(id = 'integrate', title = 'Integration method. If =1 peak limits are found through descent on the mexicanhat filtered data, if =2 the descent is done on the real data. Method 2 is very accurate but prone to noise, while method 1 is more robust to noise but less exact.', placement = 'left', trigger = 'hover'),
                          numericInput(inputId = 'max_features', label = 'maxFeatures', value = 100, min = 0),
                          bsTooltip(id = 'max_features', title = 'Maximum number of features to be defined in one bin.', placement = 'left', trigger = 'hover')
                      )
                  ),
                  box(width = 12,
                      title = 'Upload or save parameters',
                      downloadButton(outputId = 'save_params', label = 'Save parameters'),
                      fileInput(inputId = 'upload_params', label = 'Upload parameters:', multiple = F, accept = '.RData')
                  ),
                ),
                column(
                  width = 8,
                  box(
                    width = 12,
                    title = 'Results', 
                    textOutput(outputId = 'peakamount'),
                    textOutput(outputId = 'matches'),
                    plotlyOutput(outputId = 'foundpeaks')
                  ),
                  box(
                    width = 12,
                    title = 'Selected peaks',
                    verbatimTextOutput(outputId = 'vsp'),
                    uiOutput(outputId = 'mztabui')
                    # tabsetPanel(id = "mztabs", 
                    #             tabPanel('1', plotOutput('mzplot')))
                  ))
                )),
      tabItem('clustering',
              fluidRow(
                box(
                  title = 'Clustertabitem'
                )
              )),
      tabItem('statistics',
              fluidRow(
                box(
                  title = 'Statistical analysis tabitem'
                )
              )),
      tabItem('msea',
              fluidRow(
                box(
                  title = 'MSEA tabitem'
                )
              ))
    )
  ),
  footer = dashboardFooter(
    left_text = 'MetabOracle 0.3',
    right_text = 'Made by Micah Mallée'
  )
)

# Back end of the webapp is created here.
server <- function(input, output, session){
  ### functions
  {
  # Define functions
  tophits <- function(similarity_scores, limit = 5, database, splashmatches) {
    indexes_per_sample <- vector(mode = 'list', length = length(similarity_scores))
    indexes_per_sample[[1]] <- lapply(seq_along(similarity_scores[[1]]), function(x){
      top5_scores <- sort(unlist(similarity_scores[[1]][[x]]), decreasing = T)[1:limit]
      top5_matches <- vector(mode = 'list', length = length(top5_scores))
      top5_matches <- lapply(top5_scores, function(y){
        matchandscore <- list(Match = '', Score = 0)
        index1 <- grep(pattern = y, x = unlist(similarity_scores[[1]][[x]]))
        if (!is.na(splashmatches[[1]][[x]][index1])) {
          matchandscore[["Match"]] <- database[[splashmatches[[1]][[x]][index1[1]]]]
          matchandscore[["Score"]] <- round((y * 100), 2)
        } else {
          matchandscore <- NULL
        }
        matchandscore
      })
    })
  }
  
  
  ################### FIXEN matches per sample wordt niet gebruikt hahahaah
  similarities <- function(msp_query, database, SPLASH_hits) {
    #' Calculate Spectrumsimilarity per SPLASH match
    # Loop through hits
    allmatches <- vector(mode = 'list', length = length(SPLASH_hits))
    allmatches <- lapply(seq_along(SPLASH_hits), function(x){
      matches_per_sample <- vector(mode = 'list', length = length(SPLASH_hits[[x]]))
      lapply(seq_along(SPLASH_hits[[x]]), function(y){
        matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(SPLASH_hits[[x]][[y]]))
        matches_per_pseudospectrum <- sapply(SPLASH_hits[[x]][[y]], function(z){
          SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
        })
      })
    })
    # replace NaN values with 0
    allmatches <- lapply(allmatches, function(x){
      x <- lapply(x, function(y){
        y <- sapply(y, function(z){
          z <- ifelse(is.nan(z), 0, z)
        })
      })
    })
    return(allmatches)
  }
  
  matchsecondblocks <- function(querysecondblocks, databasesecondblocks) {
    #' This function takes the vectors containing the second SPLASH blocks for all pseudospectra.
    #' Matches second block with second blokcs in database per sample. Returns list with matches per pseudospectra
    matchindexes <- vector(mode = 'list', length = length(querysecondblocks))
    matchindexes <- lapply(querysecondblocks, function(x) {
      matchindexes1 <- vector(mode = 'list', length = length(x))
      matchindexes1 <- sapply(x, function(y){
        grep(pattern = y, x = databasesecondblocks)
      })
    })
  }
  
  
  getsecondblocks <- function(splashscores){
    #' This function takes all second blocks from the SPLASH codes per pseudospectrum.
    if (class(splashscores) == "character") {
      blocks <- vector(mode = 'character', length = length(splashscores))
      blocks <- sapply(splashscores, function(x){
        str_split(string = x, pattern = '-', simplify = F)[[1]][[2]]
      })
    } else {
      blocks <- vector(mode = "list", length = length(splashscores))
      blocks <- lapply(seq_along(splashscores), function(x) {
        blocks_psample <- vector(mode = 'character', length = length(splashscores[[x]]))
        blocks_psample <- sapply(splashscores[[x]], function(y) {
          str_split(string = y, pattern = '-', simplify = F)[[1]][[2]]
        })
      })
    }
    return(blocks)
  }
  
  getsplashscores <- function(msp_object) {
    #' This function calculates all SPLASH codes per pseudospectrum
    if (class(msp_object[[1]][[1]]) == "character") {
      splashscoresquery <- vector(mode = "character", length = length(msp_object))
      splashscoresquery <- sapply(seq_along(msp_object), function(y) {
        getSplash(msp_object[[y]]$pspectrum)
      })
    } else {
      splashscoresquery <- vector(mode = "list", length = length(msp_object))
      splashscoresquery <- lapply(seq_along(msp_object), function(x) {
        lapply(1:length(msp_object[[x]]), function(y) {
          getSplash(msp_object[[x]][[y]][, 1:2])
        })
      })
    }
    return(splashscoresquery)
  }
  
  split_annotate <- function(mSet) {
    #' This function splits the mSet object based on sample. 
    #' Subsequently this function groups the peaks into pseudospectra.
    f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
    f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
      gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
    })
    splitxcms <- xcms:::split.xcmsSet(smSet$xcmsSet, f = factor(f))
    annotatedxcmslist <- vector(mode = 'list', length = length(f))
    annotatedxcmslist <- lapply(splitxcms, function(x){
      x <- xsAnnotate(x)
      x <- groupFWHM(x)
    })
  }
  
  create_xchr <- function(mSet) {
    chr <- chromatogram(mSet[['onDiskData']])
    xchr <- as(chr, 'XChromatograms')
    chrompks <- mSet[["msFeatureData"]][["chromPeaks"]]
    chrompkd <- mSet[["msFeatureData"]][["chromPeakData"]]
    rt <- mSet[["msFeatureData"]][["adjustedRT"]]
    samples <- factor(chrompks[, "sample"], levels = 1:length(fileNames(mSet$onDiskData)))
    chrompks <- split.data.frame(chrompks, samples)
    chrompkd <- split.data.frame(chrompkd, samples)
    if (length(xchr) > 1) {
      for (i in 1:length(xchr)) {
        xchr[[i]]@rtime <- rt[[i]]
        xchr[[i]]@chromPeaks <- chrompks[[i]]
        xchr[[i]]@chromPeakData <- chrompkd[[i]]
      }
    } else {
      for (i in 1:length(xchr)) {
        xchr[[i]]@chromPeaks <- chrompks[[i]]
        xchr[[i]]@chromPeakData <- chrompkd[[i]]
      }
    }
    return(xchr)
  }
}
  # Create a reactive value object
  rvalues <- reactiveValues()


  
  # Dynamically set parameters
  rvalues$parameters <- reactive({
    # req(input$upload_params)
    param_file <- input$upload_params 
    if (is.null(param_file)) {
      param_initial <- SetPeakParam('general', RT_method = 'loess')
    } else {
      load(input$upload_params$datapath)
    }
    param_initial
  })
  

  # Dynamically update parameters back- and front-end
  observe({
    rvalues$param_initial <- rvalues$parameters()
    updateNumericInput(session = session, inputId = 'ppm', label = 'ppm', value = rvalues$param_initial$ppm, min = 0)
    updateNumericInput(session = session, inputId = 'noise', label = 'Noise', value = rvalues$param_initial$noise, min = 0)
    updateNumericInput(session = session, inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = rvalues$param_initial$min_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'mz_diff', label = 'mz diff', value = rvalues$param_initial$mzdiff, min = -0.01)
    updateNumericInput(session = session, inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = rvalues$param_initial$max_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'snthresh', label = 'Signal to noise threshold', value = rvalues$param_initial$snthresh, min = 0)
    updateNumericInput(session = session, inputId = 'prefilter', label = 'Prefilter', value = rvalues$param_initial$prefilter, min = 0)
    updateNumericInput(session = session, inputId = 'v_prefilter', label = 'Value of prefilter', value = rvalues$param_initial$value_of_prefilter, min = 0)
    updateSelectInput(session, inputId = 'rtmethod', label = 'Method', choices = list('loess' = 'loess', 'obiwarp' = 'obiwarp'),selected = rvalues$param_initial$RT_method)
    updateNumericInput(session, inputId = 'bw', label = 'Bandwith', value = rvalues$param_initial$bw, min = 1)
    updateNumericInput(session, inputId = 'min_fraction', label = 'minFraction', value = rvalues$param_initial$minFraction, min = 0)
    updateNumericInput(session, inputId = 'min_samples', label = 'minSamples', value = rvalues$param_initial$minSamples, min = 1)
    updateSelectInput(session, inputId = 'fitgauss', label = 'Fitgauss', choices = list('False' = 'FALSE', 'True' = 'TRUE'), selected = rvalues$param_initial$fitgauss)
    updateSelectInput(session, inputId = 'verbose_columns', label = 'Verbose columns', choices = list('False' = 'FALSE', 'True' = 'TRUE'), selected = rvalues$param_initial$verbose.columns)
    updateNumericInput(session, inputId = 'integrate', label = 'Integrate', value = rvalues$param_initial$integrate, min = 0, max = 1)
    updateNumericInput(session, inputId = 'span', label = 'Span', value = rvalues$param_initial$span, min = 0)
    updateNumericInput(session, inputId = 'max_features', label = 'maxFeatures', value = rvalues$param_initial$maxFeatures, min = 0)
    updateNumericInput(session, inputId = 'extra', label = 'Extra', value = rvalues$param_initial$extra, min = 0)
    updateNumericInput(session, inputId = 'prof_step', label = 'profStep', value = rvalues$param_initial$prof_step, min = 0)
  })
  
  
  output$file <- renderTable(input$data_input[, 1:2])
  
  # Run button (upload vs inspect) check radiobutton, then run accordingly
  observeEvent(input$run, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating plot", value = 10)
    if(input$inspect_trim == 1){
      output$inspect_plot <- renderPlot(PerformDataInspect(input$data_input$datapath))
      }
    else {
      # MSConvert_CMD <- paste0("docker run --rm -v `pwd`:`pwd` -w `pwd` chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ", input$data_input$datapath, " --mzXML")
      # system(MSConvert_CMD)
      # system('ls')
      rvalues$raw_data <- readMSData(files = input$data_input$datapath, mode = 'onDisk')
      output$inspect_plot <- renderPlot(plot(chromatogram(rvalues$raw_data)))
      }
    })
  
  # Check amount of samples. If less than 2, do not show alignment parameters
  observe({
    if (length(input$data_input$datapath) >= 0) {
      shinyjs::show(id = 'align_param_box')
    } else {
      shinyjs::hide(id = 'align_param_box')
      
    }
    
    if (input$rtmethod == 'loess') {
      shinyjs::show(id = 'loessparams')
      shinyjs::hide(id = 'obiwarpparams')
    } else if (input$rtmethod == 'obiwarp') {
      shinyjs::hide(id = 'loessparams')
      shinyjs::show(id = 'obiwarpparams')
    }
  })
  
  # Run peakdetection
  observeEvent(input$peakdetectrun, {
    params <- SetPeakParam(platform = 'general', Peak_method = 'centWave', mzdiff = input$mz_diff, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                           snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
    if (length(input$data_input$datapath) == 1) {
      mSet <- PerformPeakPicking(rvalues$raw_data, updateRawSpectraParam(params))  
    } else if(length(input$data_input$datapath) == 0) {
      showModal(modalDialog(
        title = HTML('<span style="color:#8C9398; font-size: 20px; font-weight:bold; font-family:sans-serif "> Not so fast! <span>'),
        'Please upload data first.',
        easyClose = T
      ))
      return()
    } else {
      mSet <- PerformPeakPicking(rvalues$raw_data, updateRawSpectraParam(params))
      mSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- mSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
      mSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
      mSet <- PerformPeakAlignment(mSet, param = updateRawSpectraParam(params))
      mSet <- PerformPeakFiling(mSet, param = updateRawSpectraParam(params))
      # Split mSet to nested xcmslist for further annotation
      annotatedxcmslist <- split_annotate(mSet)
      
      # Convert to MSP format
      # to.msp(object = test_data_xsannotate, file = NULL, settings = NULL, ndigit = 3, minfeat = 1, minintens = 0, intensity = "maxo", secs2mins = F)
      mspxcmslist <- lapply(X = annotatedxcmslist, to.msp, file = NULL, settings = NULL, ndigit = 3, minfeat = 1, minintens = 0, intensity = "maxo", secs2mins = F)
      querysplash <- getsplashscores(mspxcmslist)
      query_secondblocks <- getsecondblocks(querysplash)
      
      # Match secondblocks
      splashmatches <- matchsecondblocks(querysecondblocks = query_secondblocks, databasesecondblocks = mona_secondblocks)
      
      # Get similarity scores
      similarity_scores <- similarities(msp_query = mspxcmslist, database = mona_msp, SPLASH_hits = splashmatches)
      
      # Get top x matches
      bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = splashmatches)
    }
    oke <- lapply(seq_along(bestmatches), function(x){
      try(expr = {
        print(paste0(x, '  Best match: ', bestmatches[[x]][[1]][[1]][["Name"]], ' Score: ', bestmatches[[x]][[1]][["Score"]]))
      }, silent = T)
    })
    output$matches <- renderText(oke)
    rvalues$mSet <- mSet
    output$peakamount <- renderText(paste0('Amount of found peaks: ', mSet[["msFeatureData"]][["chromPeakData"]]@nrows))
    # plot plotly chromatogram with found peaks. 
    xchr <- create_xchr(mSet)
    sd <- SharedData$new(xchr)
    p <- plot_ly(x =  sd$origData()[[1]]@rtime, y = sd$origData()[[1]]@intensity, type = 'scatter', mode = 'lines', name = 'intensities', source = 'peakplot')
    for (i in 1:length(sd$origData())) {
      p <- p %>% add_trace(x =  sd$origData()[[i]]@chromPeaks[,4], y = sd$origData()[[i]]@chromPeaks[,9], 
                           type = 'bar', name = paste0('Sample ', i), text = sd$origData()[[i]]@chromPeaks[,1], 
                           hoverinfo = 'text') %>% event_register("plotly_selecting")
    }
    output$foundpeaks <- renderPlotly(p)
  })
  
  output$vsp <- renderPrint({
    d <- event_data(event = "plotly_click", priority = "event", source = 'peakplot')
    if (is.null(d)) "Click events appear here (double-click to clear)" 
    else {
      rvalues$mSet$msFeatureData$chromPeaks[which(rvalues$mSet$msFeatureData$chromPeaks[, 'rt'] == d$x),]
    }
  })
  
  output$mztabui <- renderUI({
    d <- event_data(event = "plotly_click", priority = "event", source = 'peakplot')
    if (is.null(d)) {
      NULL
    }
    else {
      pkinfo <- rvalues$mSet$msFeatureData$chromPeaks[rvalues$mSet$msFeatureData$chromPeaks[, 'rt'] == d$x,, drop = F]
      nTabs = nrow(pkinfo)
      myTabs <- lapply(paste('Tab', 1:nTabs), tabPanel, ... = plotOutput(outputId = paste('Plot', 1:nTabs)))
      do.call(tabsetPanel, myTabs)
      for (i in myTabs) {
        NULL
      }
      }
  })
  
  
  # output$mztabui <- renderUI({
  #   tabsetPanel(id = "mztabs",
  #               tabPanel('1', plotOutput('mzplot')))
  #   d <- event_data(event = "plotly_click", priority = "event", source = 'peakplot')
  #   if (is.null(d)) {
  #     NULL
  #   }
  #   else {
  #     pkinfo <- rvalues$mSet$msFeatureData$chromPeaks[rvalues$mSet$msFeatureData$chromPeaks[, 'rt'] == d$x,, drop = F]
  #     for (i in 1:ncol(pkinfo)) {
  #       insertTab(inputId = 'mztabs', target = "1", tab = tabPanel(paste0(i + 1), plotOutput(outputId = paste0("mzplot", i))))
  #       # plot(filterRt(rvalues$mSet[["onDiskData"]], rt = c(pkinfo[i, 4] - 0.001, pkinfo[i, 4] + 0.001)))
  #       # renderPlot(filterRt(rvalues$mSet[["onDiskData"]], rt = c(500, 501)))
  #     }
  #   }
  # })
  
  
  # output$mzplot <- renderPlot({
  #   d <- event_data(event = "plotly_click", priority = "event", source = 'peakplot')
  #   if (is.null(d)) {
  #     NULL
  #   }
  #   else {
  #     renderUI()
  #     pkinfo <- as.data.frame(rvalues$mSet$msFeatureData$chromPeaks[which(rvalues$mSet$msFeatureData$chromPeaks[, 'rt'] == d$x),])
  #     for (i in 1:ncol(pkinfo)) {
  #       insertTab(inputId = 'mztabs', target = "1", tab = tabPanel(paste0(i + 1), plotOutput(outputId = paste0("mzplot", i))))
  #       plot(filterRt(rvalues$mSet[["onDiskData"]], rt = c(pkinfo[4, i] - 0.001, pkinfo[4, i] + 0.001)))
  #     }
  #   }
  # })
  
  # Run automatic parameter detection and update page with new values. NOT DONE
  observeEvent(input$paramdetectrun, {
    param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter, mzdiff = input$mz_diff)
    rvalues$optimized_params <- PerformParamsOptimization(raw_data = subset_raw_data, param = param_initial, ncore = 8)
    updateNumericInput(session = session, inputId = 'ppm', label = 'ppm', value = rvalues$optimized_params$best_parameters$ppm , min = 0)
    updateNumericInput(session = session, inputId = 'noise', label = 'Noise', value = rvalues$optimized_params$best_parameters$noise, min = 0)
    updateNumericInput(session = session, inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = rvalues$optimized_params$best_parameters$min_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'mz_diff', label = 'mz diff', value = rvalues$optimized_params$best_parameters$mzdiff, min = -0.01)
    updateNumericInput(session = session, inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = rvalues$optimized_params$best_parameters$max_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'snthresh', label = 'Signal to noise threshold', value = rvalues$optimized_params$best_parameters$snthresh, min = 0)
    updateNumericInput(session = session, inputId = 'prefilter', label = 'Prefilter', value = rvalues$optimized_params$best_parameters$prefilter, min = 0)
    updateNumericInput(session = session, inputId = 'v_prefilter', label = 'Value of prefilter', value = rvalues$optimized_params$best_parameters$value_of_prefilter, min = 0)
  })
  
  # Allow upload and saving of current parameters for future use
  output$save_params <- downloadHandler(
    filename = function() {
      paste('params_', Sys.Date(), '.RData', sep = '')
    },
    content = function(file) {
      param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = input$rtmethod, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                            snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter, mzdiff = input$mz_diff)
      # save(param_initial, file = file)
      save(rvalues$mSet, file = file)
    }
  )
}







# Run app
shinyApp(ui, server)

