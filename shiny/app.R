library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(crosstalk)
library(plotly)
library(magrittr)
library(stringr)
library(MetaboAnalystR)
library(xcms)
library(metaMS)
library(OrgMassSpecR)
library(splashR)

# Increase max upload size to 400 MB
options(shiny.maxRequestSize=100*4096^2)

# Preload MoNA_DB and SPLASH hashes
mona_msp <- readRDS(file = 'data/mona_msp')
mona_splashes <- readRDS(file = 'data/mona_splashes')
# msetje <- readRDS(file = 'data/mset')
# matchmatrix <- readRDS(file = 'data/matchestable')


# Define UI.
# The front end of the webapp is created here.
ui <- dashboardPage(
  skin = 'midnight',
  header = dashboardHeader(
    title = tagList(
      tags$span(
        class = "logo-lg", "MetabOracle"
      )
    )
  ),
  # sidenavbar
  sidebar = dashboardSidebar(
    collapsed = F, 
    width = 330,
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
                         plotlyOutput(outputId = 'inspect_plot') # %>% withSpinner()
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
                  box(id = 'picking_param_box', width = 6, collapsible = T, collapsed = F,
                      style = "font-size:11px;",
                      title = 'Input parameters for peak picking',
                      numericInput(inputId = 'ppm', label = 'ppm', value = NA, min = 0),
                      bsTooltip(id = 'ppm', title = 'Maximum tolerated fluctuation of m/z value (ppm) from scan to scan - depends on the mass spectrometer accuracy', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'noise', label = 'Noise', value = 1000, min = 0),
                      bsTooltip(id = 'noise', title = 'Each centroid must be greater than the "noise" intensity value', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = 5, min = 0),
                      bsTooltip(id = 'min_peakwidth', title = 'Minimum peak width in seconds', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'mz_diff', label = 'mz diff', value = 0.05, min = -0.01),
                      bsTooltip(id = 'mz_diff', title = 'Minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = 20, min = 0),
                      bsTooltip(id = 'max_peakwidth', title = 'Maximum peak width in seconds', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'snthresh', label = 'Signal to noise threshold', value = 10, min = 0),
                      bsTooltip(id = 'snthresh', title = 'Signal to noise ratio cut-off (intensity)', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'prefilter', label = 'Prefilter', value = 3, min = 0),
                      bsTooltip(id = 'prefilter', title = 'A peak must be present in x scans with an intensity greater than the value of the prefilter', placement = 'right', trigger = 'hover'),
                      numericInput(inputId = 'v_prefilter', label = 'Value of prefilter', value = 100, min = 0),
                      bsTooltip(id = 'v_prefilter', title = 'Value of prefilter', placement = 'right', trigger = 'hover'),
                      actionButton(inputId = 'peakdetectrun', label = 'Perform peak detection', width = '100%', style = 'margin-bottom:8px;'),
                      actionButton(inputId = 'peakannotationrun', label = 'Perform peak annotation', width = '100%', style = 'margin-bottom:8px;'),
                      actionButton(inputId = 'paramdetectrun', HTML("Automatic parameter <br/>optimization"), width = '100%')
                  ),
                    box(id =  'align_param_box', width = 6, collapsible = T, collapsed = T, closable = T,
                        style = "font-size:11px;",
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
                    ),
                  box(
                    width = 12,
                    title = 'Annotation parameters', 
                    collapsible = T, 
                    collapsed = T,
                    numericInput(inputId = 'perfwhm', label = 'Full width half maximum', value = 0.6, min = 0, max = 1, step = 0.1),
                    numericInput(inputId = 'ndigit', label = 'ndigit', value = 3, min = 0, max = 100, step = 1),
                    numericInput(inputId = 'minfeat', label = 'Minimal features per pseudospectrum', value = 5, min = 0, max = 100, step = 1),
                    numericInput(inputId = 'minintens', label = 'Minimal intensity', value = 0, min = 0, max = 100000, step = 1)
                  ),
                  box(width = 12, collapsible = T, collapsed = T,
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
                  ),
                  box(
                    id = 'compoundbox',
                    collapsible = T,
                    collapsed = T, closable = T,
                    width = 12,
                    title = 'Found compounds ',
                    style = "height:500px; overflow-y: scroll; overflow-x:scroll;font-size:11px;background-color: white;",
                    class = 'cell-border stripe',
                    DTOutput(outputId = 'foundcompoundstable') %>% withSpinner()
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
    left = 'MetabOracle 0.4',
    right = 'Made by Micah Mallée'
  )
)

# Back end of the webapp is created here.
server <- function(input, output, session){
  ### functions
  {
  plot_chrom_tic_bpc <- function(raw_data) {
    plotData <- data.frame(scantime = rtime(raw_data), tic = tic(raw_data), bpc = bpi(raw_data))
    p <- plot_ly(source = "p") %>% 
      add_trace(data = plotData, x = ~scantime, y = ~tic, type = "scatter", mode = "lines", line = list(color = "rgba(0, 0, 0,0.7)", width = 0.8), 
                text = ~paste(scantime, "s"), name = "<b>Total ion chromatogram</b>") %>% 
      add_trace(data = plotData, x = ~scantime, y = ~bpc, type = "scatter", mode = "lines", line = list(color = "rgba(0, 215, 167,1)", width = 1.1), 
                text = ~paste(round(bpc, 4), "m/z"), name = "<b>Base peak chromatogram</b>") %>% 
      layout(legend = list(x = 0.7, y = 0.99), 
             xaxis = list(title = "Scan time (s)", range = c(0, max(plotData$scantime)), showspikes = TRUE, spikemode = "toaxis+across", spikesnap = "data", 
                          showline = FALSE, zeroline = FALSE, spikedash = "solid", showgrid = TRUE), 
             yaxis = list(title = "Counts", showgrid = FALSE, showticklabels = TRUE, zeroline = FALSE, showline = FALSE), hovermode = "x", showlegend = TRUE) %>% 
      event_register("plotly_click")
    return(p)
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
  
  tophits <- function(similarity_scores, limit = 5, database, splashmatches, score_cutoff = 0.8) {
    lapply(seq_along(similarity_scores), function(z){
      indexes_per_sample <- vector(mode = 'list', length = length(similarity_scores))
      indexes_per_sample <- lapply(seq_along(similarity_scores[[z]]), function(x){
        top5_scores <- sort(unlist(similarity_scores[[z]][[x]]), decreasing = T)[1:limit]
        top5_matches <- vector(mode = 'list', length = length(top5_scores))
        top5_matches <- lapply(top5_scores, function(y){
          index1 <- grep(pattern = y, x = unlist(similarity_scores[[z]][[x]]))
          if (!is.na(splashmatches[[z]][[x]][index1]) & as.numeric(y) > score_cutoff) {
            paste0('Compound: ', as.character(database[[splashmatches[[z]][[x]][index1[1]]]]["Name"]), '   Score: ', round((y * 100), 2))
          } else {
            NULL
          }
        })
      })
      names(indexes_per_sample) <- c(seq_along(indexes_per_sample))
      indexes_per_sample <- Filter(function(k) length(k) > 0, indexes_per_sample)
    })
  }
  
  similarities_thirdblocks <- function(nine_matches, msp_query, database) {
    #' Calculate Spectrumsimilarity per SPLASH match
    # Loop through hits
    if (length(nine_matches) == 1) {
      msp_query <- list(msp_query)
    }
    if (.Platform$OS.type == 'windows') {
      num_cores <- detectCores()
      cl <- makeCluster(num_cores)
      clusterExport(cl=cl, 'mona_msp')
      allmatches <- vector(mode = 'list', length = length(nine_matches))
      allmatches <- lapply(seq_along(nine_matches), function(x) {
        matches_per_sample <- vector(mode = 'list', length = length(nine_matches[[x]]))
        matches_per_sample <- parLapply(cl, seq_along(nine_matches[[x]]), function(y){
          matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(nine_matches[[x]][[y]]))
          matches_per_pseudospectrum <- sapply(nine_matches[[x]][[y]], function(z){
            OrgMassSpecR::SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
          })
        })
      })
      stopCluster(cl)
    } else {
      num_cores <- detectCores()
      allmatches <- vector(mode = 'list', length = length(nine_matches))
      allmatches <- lapply(seq_along(nine_matches), function(x) {
        matches_per_sample <- vector(mode = 'list', length = length(nine_matches[[x]]))
        matches_per_sample <- mclapply(mc.cores = num_cores, seq_along(nine_matches[[x]]), function(y){
          matches_per_pseudospectrum <- vector(mode = 'numeric', length = length(nine_matches[[x]][[y]]))
          matches_per_pseudospectrum <- sapply(nine_matches[[x]][[y]], function(z){
            OrgMassSpecR::SpectrumSimilarity(spec.top = msp_query[[x]][[y]][, 1:2], spec.bottom = database[[z]]$pspectrum, print.alignment = F, print.graphic = F)
          })
        })
      })
    }
    return(allmatches)
  }
  
  match_nines <- function(query_blocks, database_blocks) {
    #' This function takes the vectors containing the third SPLASH blocks for all pseudospectra.
    #' Matches the position of the character 9 in the the third blocks in database per sample. 
    #' Returns list with matches per pseudospectra.
    db_nines <- str_locate(database_blocks, '9')
    db_nines <- db_nines[, 1]
    query_nines <- str_locate(query_blocks, '9')
    query_nines <- query_nines[, 1]
    indexes_p_sample <- vector(mode = 'list', length = length(query_nines))
    indexes_p_sample <- lapply(seq_along(query_nines), function(y) {
      indexes <- sapply(seq_along(db_nines), function(z) ifelse(test = query_nines[[y]] == db_nines[[z]], yes = z, no = 0))
      indexes <- indexes[indexes != 0]
    })
    return(indexes_p_sample)
  }
  
  
  get_blocks <- function(splashscores, blocknr = 3){
    #' This function takes all n blocks from the SPLASH codes per pseudospectrum.
    if (class(splashscores) == "character") {
      blocks <- str_split(string = splashscores, pattern = '-')
      blocks <- sapply(blocks, function(x) x[[blocknr]])
    } else {
      blocks <- vector(mode = "list", length = length(splashscores))
      blocks <- lapply(seq_along(splashscores), function(x) {
        blocks_psample <- vector(mode = 'character', length = length(splashscores[[x]]))
        blocks_psample <- sapply(splashscores[[x]], function(y) {
          str_split(string = y, pattern = '-', simplify = F)[[1]][[blocknr]]
        })
      })
    }
    return(blocks)
  }
  
  
  get_splashscores <- function(msp_list) {
    #' This function calculates all SPLASH codes per pseudospectrum
    totalsplashscores <- vector(mode = 'list', length = length(msp_list))
    totalsplashscores <- lapply(msp_list, function(x) {
      splashscoresquery <- vector(mode = "character", length = length(x))
      splashscoresquery <- sapply(x, function(y) getSplash(y[, 1:2]))
    })
    return(totalsplashscores)
  }
  
  
  mSet2xcmsSet <-  function(mSet) {
    xs <- new("xcmsSet")
    xs@peaks <- mSet[["msFeatureData"]][["chromPeaks"]]
    rts <- list()
    if (class(mSet[[2]]) == "OnDiskMSnExp") {
      format <- "onDiskData"
    }
    else {
      format <- "inMemoryData"
    }
    rts$raw <- rtime(mSet[[format]])
    rts$corrected <- mSet[["msFeatureData"]][["adjustedRT"]]
    xs@rt <- rts
    xs@phenoData <- pData(mSet[[format]])
    xs@filepaths <- fileNames(mSet[[format]])
    xs@mslevel <- 1
    xs@scanrange <- range(scanIndex(mSet[[format]]))
    if (any(mSet$msFeatureData$chromPeakData$is_filled)) {
      fld <- which(mSet$msFeatureData$chromPeakData$is_filled)
      xs@filled <- as.integer(fld)
    }
    return(xs)
  }
  
  split_annotate <- function(mSet) {
    #' This function splits the mSet object based on sample. 
    #' Subsequently this function groups the peaks into pseudospectra.
    f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
    f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
      gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
    })
    splitxcms <- xcms:::split.xcmsSet(mSet$xcmsSet, f = factor(f))
    annotatedxcmslist <- vector(mode = 'list', length = length(f))
    annotatedxcmslist <- lapply(splitxcms, function(x){
      x <- xsAnnotate(x)
      x <- groupFWHM(x, perfwhm = 3)
    })
  }
  
  
  }
  
  # Create a reactive value object
  rvalues <- reactiveValues()
  
  # rvalues$mSet <- msetje
  
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
      # MSConvert_CMD <- paste0("docker run chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ", input$data_input$datapath, " --mzXML")
      # print(system('pwd'))
      # print(system('${PWD}'))
      # print(system('ls'))
      # system(MSConvert_CMD)
      # system('ls')
      rvalues$raw_data <- readMSData(files = input$data_input$datapath, mode = 'onDisk')
      output$inspect_plot <- renderPlotly(plot_chrom_tic_bpc(rvalues$raw_data))
      }
    })
  
  # Check amount of samples. If less than 2, do not show alignment parameters
  observe({
    if (length(input$data_input$datapath) < 2) {
      updateBox('align_param_box', action = 'remove')
      updateBox('picking_param_box', action = 'update', options = list(
        width = 12
      ))
    } else {
      updateBox('align_param_box', action = 'restore')
      updateBox('picking_param_box', action = 'update', options = list(
        width = 6
      ))
    }
  })
  
  observeEvent(input$peakannotationrun, {
    updateBox('compoundbox', action = 'toggle')
    mSet_xsannotate <- mSet2xcmsSet(rvalues$mSet)
    mSet_xsannotate <- xsAnnotate(mSet_xsannotate)
    mSet_xsannotate <- groupFWHM(mSet_xsannotate, perfwhm = input$perfwhm)
    mSet_msp <- to.msp(object = mSet_xsannotate, file = NULL, settings = NULL, ndigit = input$ndigit, minfeat = input$minfeat, minintens = input$minintens, intensity = "maxo", secs2mins = F)
    querySPLASH <- get_splashscores(msp_list = list(mSet_msp))
    full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
    full_mona_SPLASHES <- sapply(mona_msp, function(x){
      str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
    })
    query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
    database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
    SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
    similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = mSet_msp, database = mona_msp)
    bestmatches <- tophits(similarity_scores = similarity_scores, limit = 5, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = 0.8)
    # matchmatrix <- t(data.table::rbindlist(bestmatches[[1]]))
    output$foundcompoundstable <- renderDT({
      datatable(matchmatrix, class = 'cell-border stripe')
    })
  })
  
  # Run peak detection
  observeEvent(input$peakdetectrun, {
    params <- SetPeakParam(platform = 'general', Peak_method = 'centWave', mzdiff = input$mz_diff, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                           snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
    if (length(input$data_input$datapath) == 1) {
      mSet <- PerformPeakPicking(rvalues$raw_data, updateRawSpectraParam(params))
      mSet[["onDiskData"]]@phenoData@data[["sample_name"]] <- mSet[["onDiskData"]]@phenoData@data[["sampleNames"]]
      mSet[["onDiskData"]]@phenoData@data[["sampleNames"]] <- NULL
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
    }
    rvalues$mSet <- mSet
    output$peakamount <- renderText(paste0('Amount of found peaks: ', mSet[["msFeatureData"]][["chromPeakData"]]@nrows))
    p <- plot_chrom_tic_bpc(mSet$onDiskData)
    p <- p %>% add_trace(data = data.frame(mSet$msFeatureData$chromPeaks), 
                         x = ~rt, y = ~maxo, type = "scatter", mode = "markers", text = ~mz, 
                         name = paste(mSet$onDiskData@phenoData@data$sample_name, ' total peaks'))
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
    if (!is.null(d)) {
      pkinfo <- rvalues$mSet$msFeatureData$chromPeaks[rvalues$mSet$msFeatureData$chromPeaks[, 'rt'] == d$x,, drop = F]
      nTabs = nrow(pkinfo)
      myTabs <- lapply(paste('Tab', 1:nTabs), tabPanel, ... = plotOutput(outputId = paste('Plot', 1:nTabs)))
      do.call(tabsetPanel, myTabs)
      for (i in myTabs) {
        NULL
      }
    }
  })
  
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
