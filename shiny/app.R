library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyFiles)
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

# Increase max upload size to 400 MB per file
options(shiny.maxRequestSize=100*4096^2)

# Preload MoNA_DB and SPLASH hashes
mona_msp <- readRDS(file = 'data/mona_msp')
mona_splashes <- readRDS(file = 'data/mona_splashes')

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
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('home')),
      menuItem('Data inspection', tabName = 'datainspection', icon = icon('search')),
      menuItem('Peak detection', tabName = 'peakpicking', icon = icon('chart-bar')),
      menuItem('Compare samples', tabName = 'comparesamples', icon = icon('balance-scale'))
    )
  ),
  # Main content
  body = dashboardBody(
    tags$head(
      tags$style(HTML("
      .datatable_look {
      overflow-x: auto;
      overflow-y: auto;
      font-size: 12px;
      }
      thead {
      background-color: #F9F9F9;
      }
                      "))
    ),
    useShinyjs(),
    tabItems(
      tabItem('dashboard', 
              fluidRow(
                column(width = 3,
                  box(width = 12,
                    title = 'Load data',
                    helpText("Welcome to MetabOracle! Please upload your mzXML/mzML/netCDF data below."),
                    # sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
                    checkboxInput(
                      inputId = 'directory_flag',
                      label = 'Directory path?',
                      value = FALSE
                    ),
                    
                    conditionalPanel(
                      "input.directory_flag == 0",
                      shinyFilesButton(
                        id = "infile",
                        label = "Choose file(s)",
                        title = "Choose one or more files",
                        multiple = TRUE
                      )
                    ),
                    conditionalPanel(
                      "input.directory_flag == 1",
                      shinyDirButton(id = "indir", label = "Choose directory", title = "Choose a directory")
                    ),
                    verbatimTextOutput('filepaths'),
                    actionButton(inputId = 'upload', label = 'Upload')
                  )),
                column(width = 9,
                       box(
                         width = 12,
                         plotlyOutput(outputId = 'initial_plot')
                       )
                       ))),
      tabItem('datainspection',
              fluidRow(
                column(
                  width = 4,
                  box(id = 'inspection_param_box', width = 12, collapsible = T, style = 'font-size: 11px;',
                      title = 'Select region to inspect',
                      textOutput(outputId = 'mzrange'),
                      textOutput(outputId = 'rtrange'),
                      numericInput(inputId = 'min_mz', label = 'Starting value of m/z range', value = 72.9),
                      numericInput(inputId = 'max_mz', label = 'Ending value of m/z range', value = 73.1),
                      numericInput(inputId = 'min_rt', label = 'Starting value of rt range', value = 1840),
                      numericInput(inputId = 'max_rt', label = 'Ending value of rt range', value = 1880),
                      actionButton(inputId = 'runxicinspect', label = 'Show extracted ion plot'),
                      actionButton(inputId = 'runinspect', label = 'Show 3D spectrum'))
                ),
                column(width = 8,
                       box(id = 'xic_box', width = 12, title = 'Inspect figure', 
                           plotOutput(outputId = 'inspect_plot',height = '700px')))
              )),
      tabItem('peakpicking',
              fluidRow(
                column(
                  width = 4,
                  box(id = 'picking_param_box', width = 6, collapsible = T, collapsed = F, title = 'Input parameters for peak picking', style = "font-size:11px;",
                      selectInput(
                        inputId = 'peakpicking_flag',
                        label = 'centWave or matchedFilter algorithm?',choices = c(centWave = 'centWave', matchedFilter = 'matchedFilter')
                      ),
                      conditionalPanel(
                        "input.peakpicking_flag == 'matchedFilter'",
                        numericInput(inputId = 'fwhm', label = 'FWHM', value = 0.6, min = 0.1, max = 6),
                        bsTooltip(id = 'fwhm', title = ' The full width at half maximum of matched filtration gaussian model peak', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'sigma', label = 'Sigma', value = 12.74, min = 0),
                        bsTooltip(id = 'sigma', title = 'The standard deviation (width) of the matched filtration model peak', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'steps', label = 'Steps', value = 2, min = 0),
                        bsTooltip(id = 'steps', title = 'The number of bins to be merged before filtration', placement = 'right', trigger = 'hover'),
                        # numericInput(inputId = 'peakBinSize', label = 'Peak bin size', value = 0.05, min = -0.01),
                        # bsTooltip(id = 'peakBinSize', title = 'Minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap', placement = 'right', trigger = 'hover'),
                        # numericInput(inputId = 'max', label = 'Maximum peakwidth', value = 10, min = 0),
                        # bsTooltip(id = 'max', title = 'Maximum peak width in seconds', placement = 'right', trigger = 'hover')
                      ),
                      conditionalPanel(
                        "input.peakpicking_flag == 'centWave'",
                        numericInput(inputId = 'ppm', label = 'ppm', value = 25, min = 0),
                        bsTooltip(id = 'ppm', title = 'Maximum tolerated fluctuation of m/z value (ppm) from scan to scan - depends on the mass spectrometer accuracy', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'noise', label = 'Noise', value = 10000, min = 0),
                        bsTooltip(id = 'noise', title = 'Each centroid must be greater than the "noise" intensity value', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = 1, min = 0),
                        bsTooltip(id = 'min_peakwidth', title = 'Minimum peak width in seconds', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'mz_diff', label = 'mz diff', value = 0.05, min = -0.01),
                        bsTooltip(id = 'mz_diff', title = 'Minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = 10, min = 0),
                        bsTooltip(id = 'max_peakwidth', title = 'Maximum peak width in seconds', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'prefilter', label = 'Prefilter', value = 3, min = 0),
                        bsTooltip(id = 'prefilter', title = 'A peak must be present in x scans with an intensity greater than the value of the prefilter', placement = 'right', trigger = 'hover'),
                        numericInput(inputId = 'v_prefilter', label = 'Value of prefilter', value = 100, min = 0),
                        bsTooltip(id = 'v_prefilter', title = 'Value of prefilter', placement = 'right', trigger = 'hover')
                      ),
                      numericInput(inputId = 'snthresh', label = 'Signal to noise threshold', value = 100, min = 0),
                      bsTooltip(id = 'snthresh', title = 'Signal to noise ratio cut-off (intensity)', placement = 'right', trigger = 'hover'),
                      actionButton(inputId = 'peakdetectrun', label = 'Perform peak detection', width = '100%', style = 'margin-bottom:8px;'),
                      actionButton(inputId = 'createpspectra', label = HTML('Create pseudospectra'), width = '100%', style = 'margin-bottom:8px;'),
                      actionButton(inputId = 'peakannotationrun', label = 'Perform peak annotation', width = '100%', style = 'margin-bottom:8px;')
                      # actionButton(inputId = 'paramdetectrun', HTML("Automatic parameter <br/>optimization"), width = '100%')
                  ),
                    box(id = 'align_param_box', width = 6, collapsible = T,
                        style = "font-size:11px;",
                        title = 'Input parameters for peak alignment',
                        selectInput(inputId = 'rtmethod', label = 'Method', choices = list('loess' = 'loess', 'obiwarp' = 'obiwarp'),selected = 'loess'),
                        conditionalPanel(condition = "input.rtmethod == 'loess'",
                                         numericInput(inputId = 'extra', label = 'Extra', value = 1, min = 0),
                                         bsTooltip(id = 'extra', title = 'Number of "extra" peaks used to define reference peaks (or well-behaved peaks) for modeling time deviation. Number of Peaks > number of samples.', placement = 'left', trigger = 'hover'),
                                         numericInput(inputId = 'span', label = 'Span', value = 0.25, min = 0),
                                         bsTooltip(id = 'span', title = 'Degree of smoothing of the loess model. 0.2 to 1', placement = 'left', trigger = 'hover')),
                        conditionalPanel(condition = "input.rtmethod == 'obiwarp'",
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
                    bsTooltip(id = 'perfwhm', title = 'percentage of the width of the FWHM (0-1)', placement = 'right', trigger = 'hover'),
                    numericInput(inputId = 'minfeat', label = 'Minimal features per pseudospectrum', value = 5, min = 1, max = 100, step = 1),
                    bsTooltip(id = 'minfeat', title = 'Minimal amount of mass peaks per pseudospectrum', placement = 'right', trigger = 'hover'),
                    numericInput(inputId = 'minintens', label = 'Minimal intensity', value = 0.001, min = 0, max = 1, step = 0.01),
                    bsTooltip(id = 'minintens', title = 'Minimal fraction intensity per mass peak. (1 = same value as highest peak in pseudospectra)', placement = 'right', trigger = 'hover'),
                    numericInput(inputId = 'score_cutoff', label = 'Score cut-off', value = 0.6, min = 0, max = 1, step = 0.1, width = '100%'),
                    bsTooltip(id = 'score_cutoff', title = 'Minimum similarity value (0-1)', placement = 'right', trigger = 'hover'),
                    numericInput(inputId = 'hitamount', label = 'Amount of hits shown', value = 20, min = 1, max = 150, step = 5),
                    bsTooltip(id = 'hitamount', title = 'Amount of hits shown per pseudospectrum', placement = 'right', trigger = 'hover')
                  ),
                  box(width = 12, collapsible = T, collapsed = T,
                      title = 'Upload or save parameters',
                      downloadButton(outputId = 'save_params', label = 'Save parameters'),
                      fileInput(inputId = 'upload_params', label = 'Upload parameters:', multiple = F, accept = '.RData')
                  ),
                ),
                column(
                  width = 8,
                  tabBox(id = 'plottabbox', 
                         width = 12,
                         tabPanel(title = 'Detected peaks',
                                  textOutput(outputId = 'peakamount'),
                                  plotlyOutput(outputId = 'foundpeaks'))
                         ),
                  shinyjs::hidden(div(id = 'selectedbox',
                                      box(
                                          collapsible = T,
                                          width = 12,
                                          title = 'Selected',
                                          shinyjs::hidden(verbatimTextOutput(outputId = 'vsp')),
                                          shinyjs::hidden(plotOutput(outputId = 'mzspectrum')),
                                          shinyjs::hidden(plotOutput(outputId = 'alignedmzspectrum'))
                  ))),
                  uiOutput(outputId = 'compoundbox')
                  )
                )),
      tabItem(tabName = 'comparesamples', 
              fluidRow(column(width = 3,
                box(width = NULL, title = 'Select samples to use in heatmap creation', uiOutput(outputId = 'checkboxes')
              )),
              column(width = 9,
                box(width = NULL, title = 'Heatmap', plotOutput(outputId = 'heatmapfigure', height = 800))
              )
              ))
    )
  ),
  footer = dashboardFooter(
    left = 'MetabOracle 1.0',
    right = 'Made by Micah Mallée'
  )
)

# Back end of the webapp is created here.
server <- function(input, output, session){
  ### functions
  {
  heatmap_data <- function(sample_compound) {
    #' This function takes the bestmatches of x samples and returns a matrix 
    #' useable for plotting.
    all_names <- sapply(sample_compound, str_c)
    all_names_unique <- unique(unlist(all_names))
    compound_matrix1 <- matrix(nrow = length(sample_compound), ncol = length(all_names_unique))
    for (i in seq_along(sample_compound)) {
      for (j in seq_along(sample_compound[[i]])) {
        sample_index <- match(sample_compound[[i]][[j]], all_names_unique)
        compound_matrix1[i, sample_index] <- 1
      }
    }
    compound_matrix1[is.na(compound_matrix1)] <- 0
    colnames(compound_matrix1) <- all_names_unique
    rownames(compound_matrix1) <- names(sample_compound)
    return(compound_matrix1)
  }
    
  create_container <- function(sample_names) {
    #' This function takes the samples names and
    #' dynamically creates table containers for the amount of samples
    #'  and their respective pseudospectra.
    containerlist <- vector(mode = 'list', length = length(sample_names))
    containerlist <- lapply(seq_along(sample_names), function(x) {
      container_dt = withTags(table(
        class = 'datatable_look', 
        thead(
          tr(
            lapply(c(sample_names[[x]]), th, class = 'dt-center', colspan = 2, style = 'border-right: solid 1px; border-left: solid 1px;')),
          tr(lapply((c(rep(c('Compound', 'Score'), length(sample_names[[x]])))), th))
        )))
    })
    return(containerlist)
  }
    
  plot_alignment <- function(raw_data, tic_visibility = NULL, source = NULL) {
    #' This function takes the mSet data and creates a plotly figure showing the 
    #' retention time before and after alignment.
    if (is.null(tic_visibility)) {
      hovermode <- "x"
    } else {
      hovermode <- "closest"
    }
    samplenames <- raw_data[["onDiskData"]]@phenoData@data[["sample_name"]]
    fdataPerSample <- split.data.frame(raw_data[["onDiskData"]]@featureData@data, raw_data[["onDiskData"]]@featureData@data$fileIdx)
    RTPerSample <- split(raw_data[["xcmsSet"]]@rt[["raw"]], raw_data[["onDiskData"]]@featureData@data$fileIdx)
    CorrectedRT <- raw_data[["xcmsSet"]]@rt[["corrected"]]
    cols <- RColorBrewer::brewer.pal(n = length(fdataPerSample) * 2, 'Paired')
    p <- plot_ly(source = source)
    for (i in seq_along(fdataPerSample)) {
      colindex <- if(i%%2 == 0) c(-i+1, -i) else c(i, i+1)
      plotData <- data.frame(scantime = fdataPerSample[[i]]$retentionTime, tic = fdataPerSample[[i]]$totIonCurrent, bpc = fdataPerSample[[i]]$basePeakIntensity)
      p <- p %>% 
        add_trace(data = plotData, x = ~RTPerSample[[i]], y = ~bpc, type = "scatter", mode = "lines", line = list(color = cols[colindex[1]], width = 0.8), 
                  text = ~paste(scantime, "s"), visible = tic_visibility, name = paste("<b> Raw RT </b> ", samplenames[i]), legendgroup = 'Raw') %>% 
        add_trace(data = plotData, x = ~CorrectedRT[[i]], y = ~bpc, type = "scatter", mode = "lines", line = list(color = cols[colindex[1]], width = 0.8), 
                  text = ~paste(scantime, "s"), visible = tic_visibility, name = paste("<b> Corrected RT </b> ", samplenames[i]), legendgroup = 'Corrected')
    }
    p <- p %>% 
      layout(legend = list(x = 0.7, y = 0.99), 
             xaxis = list(title = "Scan time (s)", range = c(0, max(plotData$scantime)), showspikes = TRUE, spikemode = "toaxis+across", spikesnap = "data", 
                          showline = FALSE, zeroline = FALSE, spikedash = "solid", showgrid = TRUE), 
             yaxis = list(title = "Counts", showgrid = FALSE, showticklabels = TRUE, zeroline = FALSE, showline = FALSE), hovermode = "closest", showlegend = TRUE) %>% 
      event_register("plotly_click")
    return(p)
  }
    
  plot_align_spectra <- function(topSpectrum, botSpectrum) {
    #' This function takes two mass spectra and aligns them in a plotly figure.
    top_spectrum <- data.frame(mz = topSpectrum[, 1], intensity = topSpectrum[, 2])
    top_spectrum$normalized <- round((top_spectrum$intensity/max(top_spectrum$intensity)) * 100)
    fig <- plot_ly() %>% add_bars(
      y = -1*abs(top_spectrum$normalized),
      x = top_spectrum$mz,
      name = 'Pseudospectrum sample'
    )
    fig <- fig %>% add_bars(
      x = botSpectrum[,1],
      y = botSpectrum[,2],
      name = 'Pseudospectrum database'
    )
    return(fig)
  }  
    
  get_plotData_compounds <- function(bestmatches, plotData) {
    #' This function takes the bestmatches and plotData in order to create the plotData
    #' needed for the plotly figure of the detected compounds.
    compound_plotData <- vector('list', length = length(bestmatches))
    compound_plotData <- lapply(seq_along(bestmatches), function(x) {
      pspectra_nr <- as.integer(sub('.*Pseudospectrum_', '', names(bestmatches[[x]])))
      compound_plotData_per_sample <- plotData[[x]][as.integer(pspectra_nr),]
      top_compounds_per_sample <- vector(mode = 'character', length = length(bestmatches[[x]]))
      top_compounds_per_sample <- sapply(seq_along(bestmatches[[x]]), function(y) {
        bestmatches[[x]][[y]][[1]][1]
      })
      compound_plotData_per_sample <- cbind(compound_plotData_per_sample, Compound = top_compounds_per_sample)
      compound_plotData_per_sample
    })
  }
    
  plotdata_pseudospectra <- function(msps){
    #' This function takes the MSP data in order to create the plotData
    #' needed for the plotly figure of the created pseudospectra.
    plotData <- vector(mode = 'list', length = length(msps))
    for (i in seq_along(msps)){
      pseudospectrum <- lapply(seq_along(msps[[i]]), function(x){
        rt <- round(median(msps[[i]][[x]][, 3]), 4)
        intense <- max(msps[[i]][[x]][, 2])
        list(rt = rt, intense = intense, pspectra_id = x, sample_nr = i)
      })
      plotData[[i]] <- pseudospectrum
    }
    plotData <- lapply(plotData, data.table::rbindlist, use.names = T)
    return(plotData)
  }
    
  plot_chrom_tic_bpc <- function(raw_data, tic_visibility = NULL, source = NULL) {
    #' This function takes the raw data and creates a plotly figure showing the 
    #' BPC and TIC.
    if (is.null(tic_visibility)) {
      hovermode <- "x"
    } else {
      hovermode <- "closest"
    }
    samplenames <- raw_data@phenoData@data[["sample_name"]]
    fdataPerSample <- split.data.frame(raw_data@featureData@data, raw_data@featureData@data$fileIdx)
    cols <- RColorBrewer::brewer.pal(n = length(fdataPerSample) * 2, 'Paired')
    p <- plot_ly(source = source)
    for (i in seq_along(fdataPerSample)) {
      colindex <- if(i%%2 == 0) c(-i+1, -i) else c(i, i+1)
      plotData <- data.frame(scantime = fdataPerSample[[i]]$retentionTime, tic = fdataPerSample[[i]]$totIonCurrent, bpc = fdataPerSample[[i]]$basePeakIntensity)
      p <- p %>% 
        add_trace(data = plotData, x = ~scantime, y = ~tic, type = "scatter", mode = "lines", line = list(color = cols[colindex[1]], width = 0.8), 
                  text = ~paste(scantime, "s"), visible = tic_visibility, name = paste("<b>Total ion chromatogram</b> ", samplenames[i])) %>% 
        add_trace(data = plotData, x = ~scantime, y = ~bpc, type = "scatter", mode = "lines", line = list(color = cols[colindex[2]], width = 1.1), 
                  text = ~paste(round(bpc, 4), "m/z"), name = paste("<b>Base peak chromatogram</b> ", samplenames[i]))
    }
    p <- p %>% 
      layout(legend = list(x = 0.7, y = 0.99), 
             xaxis = list(title = "Scan time (s)", range = c(0, max(plotData$scantime)), showspikes = TRUE, spikemode = "toaxis+across", spikesnap = "data", 
                          showline = FALSE, zeroline = FALSE, spikedash = "solid", showgrid = TRUE), 
             yaxis = list(title = "Counts", showgrid = FALSE, showticklabels = TRUE, zeroline = FALSE, showline = FALSE), hovermode = "closest", showlegend = TRUE) %>% 
      event_register("plotly_click")
    return(p)
  }
  
  tophits <- function(similarity_scores, limit = 5, database, splashmatches, score_cutoff = 0.8) {
    #' This function takes the similarity scores, database, and splashmatches in order to create 
    #' a nested list containing the tophits per pseudospectrum per sample.
    totalmatches <- vector('list', length(similarity_scores))
    totalmatches <- lapply(seq_along(similarity_scores), function(z){
      indexes_per_sample <- vector(mode = 'list', length = length(similarity_scores[[z]]))
      indexes_per_sample <- lapply(seq_along(similarity_scores[[z]]), function(x){
        top5_scores <- sort(unlist(similarity_scores[[z]][[x]]), decreasing = T)[1:limit]
        top5_matches <- vector(mode = 'list', length = length(top5_scores))
        top5_matches <- lapply(top5_scores, function(y){
          index1 <- grep(pattern = y, x = unlist(similarity_scores[[z]][[x]]))
          if (!is.na(splashmatches[[z]][[x]][index1]) & as.numeric(y) > score_cutoff) {
            c(as.character(database[[splashmatches[[z]][[x]][index1[1]]]]["Name"]), round((y * 100), 2), database[[splashmatches[[z]][[x]][index1[1]]]]["pspectrum"])
          } else {
            NULL
          }
        })
      })
      names(indexes_per_sample) <- c(sprintf("Pseudospectrum_%d", seq.int(indexes_per_sample)))
      indexes_per_sample <- Filter(Negate(function(x) is.null(unlist(x))), indexes_per_sample)
    })
    return(totalmatches)
  }
  
  
  similarities_thirdblocks <- function(nine_matches, msp_query, database) {
    #' Calculate SpectrumSimilarity per SPLASH match
    # Loop through hits
    if (.Platform$OS.type == 'windows') {
      num_cores <- detectCores()
      cl <- makeCluster(num_cores)
      clusterExport(cl=cl, c('database', 'msp_query'), envir = environment())
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
    #' This function converts an mSet object into an xcmsSet object
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
  
  split_mSet <- function(mSet) {
    #' This function splits the mSet object based on sample into xcmsSets. 
    f <- vector(mode = 'character', length = length(mSet$xcmsSet@phenoData$sample_name))
    f <- sapply(mSet$xcmsSet@phenoData$sample_name, function(x) {
      gsub(x = x, pattern = " ", replacement =  ".", fixed = T)
    })
    splitxcms <- xcms:::split.xcmsSet(mSet$xcmsSet, f = factor(f))
    for(i in seq_along(splitxcms)) {splitxcms[[i]]@peaks[, 11] <- i}
    return(splitxcms)
  }
  
  annotate_xcmslist <- function(xcmslist, perfwhm = 0.6){
    #' This function groups the peaks into pseudospectra per sample.
    annotatedxcmslist <- vector(mode = 'list', length = length(xcmslist))
    annotatedxcmslist <- lapply(xcmslist, function(x){
      x <- xsAnnotate(x)
      x <- groupFWHM(x, perfwhm = perfwhm)
    })
    return(annotatedxcmslist)
  }
  }
  
  # Create a reactive value object
  rvalues <- reactiveValues()
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  shinyFileChoose(
    input,
    'infile',
    roots = volumes
  )
  
  shinyDirChoose(
    input,
    'indir',
    roots = volumes
  )
  
  
  output$filepaths <- renderText({
    #' Here the filepath of the uploaded directory or file(s) is saved.
    if (!is.integer(input$infile)) {
      filelocation <- parseFilePaths(volumes, input$infile)
      rvalues$data_input <- filelocation
      rvalues$dir_or_file <- length(filelocation$datapath)
      rvalues$data_input$datapath
    } else if (!is.integer(input$indir)) {
      dirlocation <- parseDirPath(volumes, input$indir)
      rvalues$data_input <- dirlocation
      dirlocation
    } else {
      'No files have been selected'
    }
  })
  
  
  #' Dynamically change the UI depending on the amount of samples uploaded.
  #' If less than 2 samples are uploaded, the alignment parameters are not shown.
  observeEvent(rvalues$data_input, {
    if (is.null(rvalues$dir_or_file)){
      shinyjs::show('align_param_box')
      updateBox(id = 'picking_param_box', action = 'update', options = list(width = 6))
    } else if (!is.null(rvalues$dir_or_file) & rvalues$dir_or_file == 1) {
      shinyjs::hide('align_param_box')
      updateBox(id = 'picking_param_box', action = 'update', options = list(width = 12))
    } else if(!is.null(rvalues$dir_or_file) & rvalues$dir_or_file > 1) {
      shinyjs::show('align_param_box')
      updateBox(id = 'picking_param_box', action = 'update', options = list(width = 6))
    }
  })
  
  
  # Dynamically set parameters
  #' When parameters are uploaded, this piece of code changes the shown parameters.
  #' If no parameters are uploaded, standard parameters are shown.
  rvalues$parameters <- reactive({
    req(input$upload_params)
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
    updateNumericInput(inputId = 'fwhm', label = 'FWHM', value = rvalues$param_initial$fwhm, min = 0.1, max = 6)
    updateNumericInput(inputId = 'sigma', label = 'Sigma', value = rvalues$param_initial$sigma, min = 0)
    updateNumericInput(inputId = 'steps', label = 'Steps', value = rvalues$param_initial$steps, min = 0)
    updateSelectInput(session = session, inputId = 'peakpicking_flag', selected = rvalues$param_initial$Peak_method)
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
  
  # This observeEvent handles the data upload
  observeEvent(input$upload, {
    req(rvalues$data_input)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating plot", value = 10)
    if(is.null(rvalues$dir_or_file)) {
      rvalues$raw_data <- ImportRawMSData(foldername = rvalues$data_input, mode = 'onDisk', ncores = detectCores(), plotSettings = SetPlotParam(Plot = F))
    } else {
      sample_names <- c(rvalues$data_input$datapath)
      sample_names <- gsub("^.*/", "", sample_names)
      rvalues$raw_data <- readMSData(files = c(rvalues$data_input$datapath), mode = 'onDisk')
      rvalues$raw_data@phenoData@data[["sample_name"]] <- sample_names
      rvalues$raw_data@phenoData@data[["sampleNames"]] <- NULL
    }
    output$initial_plot <- renderPlotly(plot_chrom_tic_bpc(rvalues$raw_data, source = NULL))
    })
  
  # This observeEvent runs when the 'inspect' button is clicked. Creates a 3D plot.
  observeEvent(input$runinspect, {
    req(rvalues$data_input)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating plot", value = 10)
    output$inspect_plot <- renderPlot(PerformDataInspect(rvalues$data_input$datapath[1]))
  })
  
  # Dynamically show the mz range
  output$mzrange <- renderText({
    req(rvalues$raw_data)
    mzrange <- c(min(rvalues$raw_data@featureData@data[["lowMZ"]]), max(rvalues$raw_data@featureData@data[["highMZ"]]))
    paste('MZ values range from: ', mzrange[1], ' to ', mzrange[2])
  })
  
  # Dynamically show the rt range
  output$rtrange <- renderText({
    req(rvalues$raw_data)
    rtrange <- c(min(rvalues$raw_data@featureData@data[["retentionTime"]]), max(rvalues$raw_data@featureData@data[["retentionTime"]]))
    paste('Retention time ranges from: ', rtrange[1], ' to ', rtrange[2])
  })
  
  # This observeEvent runs when the 'XICinspect' button is clicked. Creates a XIC plot.
  observeEvent(input$runxicinspect, {
    req(input$min_mz, input$max_mz, input$min_rt, input$max_rt, rvalues$raw_data)
    withProgress(message = 'Plotting...', {
      rvalues$raw_data %>%
        filterRt(rt = c(input$min_rt, input$max_rt)) %>%
        filterMz(mz = as.numeric(c(input$min_mz, input$max_mz))) -> p_xic
      output$inspect_plot <- renderPlot(plot(p_xic, type = 'XIC'))
    })
  })


  # Observes the current selected table of compounds.
  observe({
    rvalues$tableindex <- match(input$compoundbox, names(rvalues$matchmatrices))
    rvalues$currenttable <- paste0('table', rvalues$tableindex, '_cell_clicked')
  })
  
  
  #' When the current selected table of compounds is clicked on, 
  #' this shows the aligned mass spectra for the selected cell.
  observeEvent(input[[rvalues$currenttable]], {
    rank <- input[[rvalues$currenttable]]$row
    pspectra <- input[[rvalues$currenttable]]$col + 1
    pspectra_id <- colnames(rvalues$matchmatrices[[rvalues$tableindex]])[pspectra]
    if (!is.null(pspectra) & !is.null(rank)) {
      shinyjs::show(id = 'selectedbox', anim = T)
      shinyjs::show(id = 'alignedmzspectrum', anim = T, animType = 'fade', asis = T)
      shinyjs::hide(id = 'vsp')
      shinyjs::hide(id = 'mzspectrum')
      pspectra <- as.integer(str_split(pspectra_id, '_')[[1]][2])
      output$alignedmzspectrum <- renderPlot(OrgMassSpecR::SpectrumSimilarity(spec.top = rvalues$mSet_msp[[input$compoundbox]][[pspectra]],
                                                                              spec.bottom = rvalues$bestmatches_pspectra[[rvalues$tableindex]][[pspectra_id]][[rank]]$pspectrum,
                                                                              top.label = 'Detected compound', bottom.label = 'Database compound',
                                                                              print.alignment = F))
    }
  })
  
  
  # Run peak detection
  observeEvent(input$peakdetectrun, {
    withProgress(message = 'Running peak detection', {
      tryCatch({
        if (input$peakpicking_flag == 'centWave') {
          params <- SetPeakParam(platform = 'general', Peak_method = 'centWave', mzdiff = input$mz_diff, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                 snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
        } else if (input$peakpicking_flag == 'matchedFilter') {
          params <- SetPeakParam(platform = 'general', Peak_method = input$peakpicking_flag, mzdiff = input$mz_diff,
                                 snthresh = input$snthresh, fwhm = input$fwhm, sigma = input$sigma, steps = input$steps, peakBinSize = input$peakBinSize)
        }
        if(is.null(rvalues$raw_data)) {
          showModal(modalDialog(
            title = HTML('<span style="color:#8C9398; font-size: 20px; font-weight:bold; font-family:sans-serif "> Not so fast! <span>'),
            'Please upload data first.',
            easyClose = T
          ))
          return()
        }
        if (is.null(rvalues$dir_or_file)){
          mSet <- PerformPeakProfiling(rvalues$raw_data, params, plotSettings = SetPlotParam(Plot = F))
          xcmslist <-  split_mSet(mSet = mSet)
          rvalues$xcmslist <- xcmslist
        } else if (!is.null(rvalues$dir_or_file) & rvalues$dir_or_file == 1) {
          mSet <- PerformPeakPicking(rvalues$raw_data, updateRawSpectraParam(params))
          xcmsSet <- mSet2xcmsSet(mSet)
          p <- plot_chrom_tic_bpc(mSet$onDiskData, tic_visibility = 'legendonly', source = 'p')
          samplename <- list(xcmsSet@phenoData["sample_name"][[1]])
          plotData_peaks <- list(data.frame(xcmsSet@peaks))
          rvalues$xcmsSet <- xcmsSet
        } else if(!is.null(rvalues$dir_or_file) & rvalues$dir_or_file > 1) {
          mSet <- PerformPeakPicking(rvalues$raw_data, updateRawSpectraParam(params))
          mSet <- PerformPeakAlignment(mSet, param = updateRawSpectraParam(params))
          mSet <- PerformPeakFiling(mSet, param = updateRawSpectraParam(params))
          xcmslist <-  split_mSet(mSet = mSet)
          rvalues$xcmslist <- xcmslist
        }
        if (!is.null(rvalues$xcmslist)) {
          p <- plot_chrom_tic_bpc(mSet$onDiskData, tic_visibility = 'legendonly', source = 'p')
          plotData_peaks <- vector(mode = 'list', length = length(xcmslist))
          for (i in seq_along(xcmslist)) {
            plotData_peaks[[i]] <- as.data.frame(xcmslist[[i]]@peaks)
          }
          for (i in plotData_peaks) { i$rt <- round(i$rt, 4)}
          samplename <- sapply(xcmslist, function(x) x@phenoData[["sample_name"]])
          if (is.null(rvalues$checked[[2]])) {
            appendTab(inputId = 'plottabbox', tab = tabPanel(title = 'Alignment results', plotlyOutput('alignmentresults')), select = T)
            rvalues$checked[[2]] <- 2
          }
          p_aligned <- plot_alignment(raw_data = mSet, tic_visibility = NULL, source = 'p_aligned')
          output$alignmentresults <- renderPlotly(p_aligned)
        }
        rvalues$plotData_peaks <- plotData_peaks
        for (i in seq_along(plotData_peaks)) {
          p <- p %>% add_trace(data = plotData_peaks[[i]],
                               x = ~rt, y = ~maxo, type = "scatter", mode = "markers", text = ~mz,
                               name = paste(samplename[i], ' total peaks'))
        }
        mSet$onDiskData@phenoData@data$sample_name <- rvalues$raw_data@phenoData@data[["sample_name"]]
        rvalues$mSet <- mSet
        output$peakamount <- renderText(paste0('Amount of found peaks: ', mSet[["msFeatureData"]][["chromPeakData"]]@nrows))
        output$foundpeaks <- renderPlotly(p)
      }, error = function(e) {
        showModal(modalDialog(
          title = HTML('<span style="color:#8C9398; font-size: 20px; font-weight:bold; font-family:sans-serif "> No peaks found <span>'),
          'Try changing the parameters, 
          SN treshold might be too high.',
          easyClose = T
        ))
      })
      })
  })
  
  # This observeEvent shows the peakinfo when clicked on.
  observeEvent(event_data(event = "plotly_click", priority = "event", source = 'p'), {
    shinyjs::show(id = 'selectedbox', anim = T)
    d <- event_data(event = "plotly_click", priority = "event", source = 'p')
    plotData_peaks1 <- data.table::rbindlist(rvalues$plotData_peaks, use.names = T)
    if (is.null(d)) {
      return(NULL)
    } else {
      plotData_peaks1 %>% filter(round(rt, 4) %in% round(d$x, 4)) -> peakdata
      shinyjs::show(id = 'vsp', anim = T, animType = 'fade')
      shinyjs::hide(id = 'mzspectrum', anim = T, animType = 'slide')
      output$vsp <- renderPrint(peakdata)
    }
  })

  rvalues$checked <- list(NULL, NULL, NULL)
  
  # Runs the creation of pseudospectra.
  observeEvent(input$createpspectra, {
    if (is.null(rvalues$dir_or_file) | rvalues$dir_or_file > 1){
      xcmslist <-  annotate_xcmslist(xcmslist = rvalues$xcmslist, perfwhm = input$perfwhm)
      mSet_msp <- lapply(xcmslist, to.msp, file = NULL, settings = NULL, minfeat = input$minfeat, minintens = input$minintens, intensity = "maxo", secs2mins = F)
    } else if (!is.null(rvalues$dir_or_file) & rvalues$dir_or_file == 1) {
      mSet_xsannotate <- xsAnnotate(rvalues$xcmsSet)
      mSet_xsannotate <- groupFWHM(mSet_xsannotate, perfwhm = input$perfwhm)
      mSet_msp <- to.msp(object = mSet_xsannotate, file = NULL, settings = NULL, minfeat = input$minfeat, minintens = input$minintens, intensity = "maxo", secs2mins = F)
      mSet_msp <- list(mSet_msp)
      names(mSet_msp) <- rvalues$xcmsSet@phenoData$sample_name
    }
    if (is.null(rvalues$checked[[1]])) {
      appendTab(inputId = 'plottabbox', tab = tabPanel(title = 'Pseudospectra', textOutput(outputId = 'pspectra_amount'), plotlyOutput('foundpseudospectra')), select = T)
      rvalues$checked[[1]] <- 1
    } 
    rvalues$mSet_msp <- mSet_msp
    p1 <- plot_chrom_tic_bpc(rvalues$mSet$onDiskData, tic_visibility = 'legendonly', source = 'p1')
    plotData <- plotdata_pseudospectra(msps = mSet_msp)
    rvalues$plotData <- plotData
    for (i in seq_along(plotData)) {
      p1 <- p1 %>% add_trace(data = plotData[[i]], x = ~rt, y = ~intense, type = "scatter", mode = "markers", text = ~pspectra_id, name = paste('Pseudospectra: ', names(mSet_msp)[i]))
    }
    output$pspectra_amount <- renderPrint(lengths(mSet_msp))
    output$foundpseudospectra <- renderPlotly(p1)
  })

  
  # When a compound is clicked on, this piece of code shows the aligned spectra
  observeEvent(event_data(event = "plotly_click", source = "p2"), {
    shinyjs::show(id = 'selectedbox', anim = T)
    d2 <- event_data(event = "plotly_click", source = "p2")
    plotData_compounds1 <- data.table::rbindlist(rvalues$plotData_compounds, use.names = T)
    if (nrow(plotData_compounds1) != 0) {
      plotData_compounds1 %>%
        filter(rt %in% d2$x) -> ja
    } else {
      return(NULL)
    }
    if (nrow(ja) == 0) {
      shinyjs::hide(id = 'alignedmzspectrum')
      return(NULL)
    }
    shinyjs::show(id = 'alignedmzspectrum', anim = T, animType = 'fade', asis = T)
    shinyjs::hide(id = 'vsp')
    shinyjs::hide(id = 'mzspectrum')
    output$alignedmzspectrum <- renderPlot(OrgMassSpecR::SpectrumSimilarity(spec.top = rvalues$mSet_msp[[ja$sample_nr]][[ja$pspectra_id]], 
                                                                            spec.bottom = rvalues$bestmatches_pspectra[[1]][[ja$pspectra_id]][[1]]$pspectrum, 
                                                                            top.label = 'Detected compound', bottom.label = 'Database compound', 
                                                                            print.alignment = F))
  })
  
  
  
  # When a pseudospectrum is clicked on, this piece of code shows the underlying 
  # mass spectrum.
  observeEvent(event_data(event = "plotly_click", source = "p1"), {
    shinyjs::show(id = 'selectedbox', anim = T)
    d1 <- event_data(event = "plotly_click", source = "p1")
    plotData1 <- data.table::rbindlist(rvalues$plotData, use.names = T)
    if (nrow(plotData1) != 0) {
      plotData1 %>%
        filter(rt %in% d1$x) -> ja
    } else {
      return(NULL)
    }
    if (nrow(ja) == 0) {
      shinyjs::hide(id = 'mzspectrum')
      return(NULL)
      }
    shinyjs::show(id = 'mzspectrum', anim = T, animType = 'fade', asis = T)
    shinyjs::hide(id = 'vsp')
    output$mzspectrum <- renderPlot(plotPseudoSpectrum(rvalues$mSet_msp[[ja$sample_nr]][[ja$pspectra_id]][, 1:2]))
  })
  
  # Observes when the peak annotation button is clicked. Runs peak annotation and outputs plot
  observeEvent(input$peakannotationrun, {
    if (is.null(rvalues$mSet_msp)) {
      showModal(modalDialog(
        title = HTML('<span style="color:#8C9398; font-size: 20px; font-weight:bold; font-family:sans-serif "> Not so fast! <span>'),
        'Please run peak detection and pseudospectra creation first.',
        easyClose = T
      ))
      return()
    }
    withProgress(message = 'Running peak annotation', {
      full_mona_SPLASHES <- vector(mode = 'character', length = length(mona_msp))
      full_mona_SPLASHES <- sapply(mona_msp, function(x){
        str_extract(string = x$Comments, pattern = regex('SPLASH=splash10................'))
      })
      querySPLASH <- get_splashscores(msp_list = rvalues$mSet_msp)
      query_thirdblocks <- lapply(querySPLASH, get_blocks, blocknr = 3)
      database_thirdblocks <- get_blocks(splashscores = full_mona_SPLASHES, blocknr = 3)
      SPLASH_matches <- lapply(query_thirdblocks, match_nines, database_blocks = database_thirdblocks)
      similarity_scores <- similarities_thirdblocks(nine_matches = SPLASH_matches, msp_query = rvalues$mSet_msp, database = mona_msp)
      for(i in seq_along(similarity_scores)) {names(similarity_scores[[i]]) <- sprintf("Pseudospectrum_%d", seq.int(similarity_scores[[i]]))}
      bestmatches <- tophits(similarity_scores = similarity_scores, limit = input$hitamount, database = mona_msp, splashmatches = SPLASH_matches, score_cutoff = input$score_cutoff)
      names(bestmatches) <- names(rvalues$mSet[["onDiskData"]]@phenoData@data[["sample_name"]])
      bestmatches_pspectra <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[3])))
      rvalues$bestmatches_pspectra <- bestmatches_pspectra
      bestmatches <- lapply(bestmatches, function(x) lapply(x, function(y) lapply(y, function(z) z[-3])))
      rvalues$bestmatches <- bestmatches
      matchmatrices <- vector(mode = 'list', length = length(bestmatches))
      names(matchmatrices) <- rvalues$mSet[["onDiskData"]]@phenoData@data[["sample_name"]]
      sample_names <- vector(mode = 'list', length(bestmatches))
      for (i in seq_along(bestmatches)) {
        matchmatrices[[i]] <- t(data.table::rbindlist(bestmatches[[i]]))
        colnames(matchmatrices[[i]]) <- sort(rep(names(bestmatches[[i]]), times = 2))
        sample_names[[i]] <- names(bestmatches[[i]])
      }
      containerlist <- create_container(sample_names)
      
      output$compoundbox <- renderUI({
        ntabs = length(matchmatrices)
        
        myTabs = lapply(seq_len(ntabs), function(x){
          
          tabPanel(paste0(names(matchmatrices[x])), dataTableOutput(outputId = paste0('table', x)) , class = "datatable_look") 
        })
        do.call(what = shinydashboard::tabBox, args = c(myTabs, list(id = 'compoundbox', width = 12)))
      })
      
      lapply(seq_len(length(bestmatches)), function(i) {
        output[[paste0("table",i)]] <- DT::renderDataTable(
          matchmatrices[[i]], class = 'cell-border stripe', extensions = 'Buttons',
          options = list(dom = 'Bfrtip', buttons = list('copy', 'print', list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'), text = 'Download')),
                         paging = F, searching = F),
          filter = list(position = 'top', clear = F), rownames = F,
          container = containerlist[[i]], selection = 'single'
        )
      })  
      rvalues$matchmatrices <- matchmatrices
      plotData_compounds <- get_plotData_compounds(bestmatches = bestmatches, plotData = plotdata_pseudospectra(rvalues$mSet_msp))
      rvalues$plotData_compounds <- plotData_compounds
      p2 <- plot_chrom_tic_bpc(rvalues$mSet$onDiskData, tic_visibility = 'legendonly', source = 'p2')
      for (i in seq_along(plotData_compounds)) {
        p2 <- p2 %>% add_trace(data = plotData_compounds[[i]], x = ~rt, y = ~intense, type = "scatter", mode = "markers", text = ~Compound, name = paste("Compounds: ", names(rvalues$mSet_msp)[i]))
      }
      if (is.null(rvalues$checked[[3]])) {
        appendTab(inputId = 'plottabbox', tab = tabPanel(title = 'Compounds detected', plotlyOutput('foundcompounds')), select = T)
        rvalues$checked[[3]] <- 3
      }
      saveRDS(similarity_scores, 'sim_scores')
      saveRDS(rvalues$mSet, 'mSet')
      output$foundcompounds <- renderPlotly(p2)
    })
  })
  
  
  # Run automatic parameter detection and update page with new values.
  observeEvent(input$paramdetectrun, {
    if(!is.null(rvalues$dir_or_file) | is.null(rvalues$data_input)) {
      showModal(modalDialog(
        title = HTML('<span style="color:#8C9398; font-size: 20px; font-weight:bold; font-family:sans-serif "> Not so fast! <span>'),
        'Please upload data first.',
        easyClose = T
      ))
      return()
    }
    withProgress({
      memory_raw_data <- ImportRawMSData(foldername = rvalues$data_input, mode = 'inMemory', ncores = detectCores(), plotSettings = SetPlotParam(Plot = F))
      param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                    snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter, mzdiff = input$mz_diff)
      rvalues$optimized_params <- PerformParamsOptimization(raw_data = memory_raw_data, param = param_initial, ncore = detectCores())
      updateNumericInput(session = session, inputId = 'ppm', label = 'ppm', value = rvalues$optimized_params$best_parameters$ppm , min = 0)
      updateNumericInput(session = session, inputId = 'noise', label = 'Noise', value = rvalues$optimized_params$best_parameters$noise, min = 0)
      updateNumericInput(session = session, inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = rvalues$optimized_params$best_parameters$min_peakwidth, min = 0)
      updateNumericInput(session = session, inputId = 'mz_diff', label = 'mz diff', value = rvalues$optimized_params$best_parameters$mzdiff, min = -0.01)
      updateNumericInput(session = session, inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = rvalues$optimized_params$best_parameters$max_peakwidth, min = 0)
      updateNumericInput(session = session, inputId = 'snthresh', label = 'Signal to noise threshold', value = rvalues$optimized_params$best_parameters$snthresh, min = 0)
      updateNumericInput(session = session, inputId = 'prefilter', label = 'Prefilter', value = rvalues$optimized_params$best_parameters$prefilter, min = 0)
      updateNumericInput(session = session, inputId = 'v_prefilter', label = 'Value of prefilter', value = rvalues$optimized_params$best_parameters$value_of_prefilter, min = 0)
      saveRDS(rvalues$optimized_params, 'optimized_params_rhino')
    })
  })
  
  # Allow upload and saving of current parameters for future use
  output$save_params <- downloadHandler(
    filename = function() {
      paste('params_', Sys.Date(), '.RData', sep = '')
    },
    content = function(file) {
      if (input$peakpicking_flag == 'centWave') {
        param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', mzdiff = input$mz_diff, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                               snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
      } else if (input$peakpicking_flag == 'matchedFilter') {
        param_initial <- SetPeakParam(platform = 'general', Peak_method = input$peakpicking_flag, mzdiff = input$mz_diff,
                               snthresh = input$snthresh, fwhm = input$fwhm, sigma = input$sigma, steps = input$steps, peakBinSize = input$peakBinSize)
      }
      save(param_initial, file = file)
    }
  )
  
  # Heatmap creation
  ## Checks if annotation has been performed. If yes, shows selectable samples.
  ## once user clicks the run button, heatmap is generated and shown.
  output$checkboxes <- renderUI({
    req(rvalues$mSet)
    outputs <- tagList()
    choice <- rvalues$mSet[["onDiskData"]]@phenoData@data[["sample_name"]]
    outputs[[1]] <- checkboxGroupInput(inputId = 'checkboxsamples', label = 'Select samples', choices = choice, selected = NULL)
    outputs[[2]] <- actionButton(inputId = 'runheatmap', label = 'Create heatmap with selected samples')
    outputs
  })
  
  rvalues$mSet <- readRDS('mSet')
  rvalues$bestmatches <- readRDS('bestmatches')
  
  observeEvent(input$runheatmap, {
    withProgress({
      bestmatches <- rvalues$bestmatches
      top_compounds <- lapply(seq_along(bestmatches), function(x){
        top_compounds_per_sample <- sapply(seq_along(bestmatches[[x]]), function(y) {
          bestmatches[[x]][[y]][[1]][1]
        })
      })
      names(top_compounds) <- rvalues$mSet[["onDiskData"]]@phenoData@data[["sample_name"]]
      selected_samples_heatmap_data <- heatmap_data(top_compounds[input$checkboxsamples])
      output$heatmapfigure <- renderPlot(pheatmap::pheatmap(mat = t(selected_samples_heatmap_data), scale = 'none', angle_col = 45, fontsize_row = 8, cluster_rows = F, 
                                                            color = RColorBrewer::brewer.pal(n = 3, name = "Set1")))
    })
  })
  
  
}


# Run app
shinyApp(ui, server)
