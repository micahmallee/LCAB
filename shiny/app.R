library(shiny)
library(MetaboAnalystR)
library(shinydashboard)
library(shinydashboardPlus)
library(xcms)
library(magrittr)
library(shinycssloaders)
library(shinyBS)
library(shinyjs)

options(shiny.maxRequestSize=100*1024^2)

# Define UI
ui <- dashboardPagePlus(
  skin = 'midnight',
  header = dashboardHeaderPlus(
    title = tagList(
      tags$span(
        class = "logo-lg", "MetabOracle"
      )
    )
  ),
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
    
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem('dashboard', 
              fluidRow(
                box(
                  title = 'Load data',
                  helpText("Welcome to MetabOracle! Please upload your mzXML/mzML/netCDF data below."),
                  fileInput(inputId = 'data_input', label = 'Select your file(s):', multiple = T, accept = c('.mzXML', 'mzxml', '.CDF', '.cdf', '.mzml', 'mzML')),
                  sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
                  radioButtons(inputId = 'inspect_trim', label = 'Inspect or upload data.',
                               choices = list('Inspect' = 1,
                                              'Upload' = 2),
                               selected = 1, inline = T),
                  actionButton(inputId = 'run', label = 'Run')
                ),
                box(
                  imageOutput(outputId = 'inspect_plot') # %>% withSpinner()
                ),
                box(
                  title = 'File information:',
                  tableOutput(outputId = 'file'),
                  textOutput('paths')
                ))),
      tabItem('peakpicking',
              fluidRow(
                box(width = 2,
                  title = 'Input parameters for peak picking',
                  numericInput(inputId = 'ppm', label = 'ppm', value = NA, min = 0),
                  bsTooltip(id = 'ppm', title = 'Maximum mass deviation', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'noise', label = 'Noise', value = 1000, min = 0),
                  bsTooltip(id = 'noise', title = 'Noise value', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = 5, min = 0),
                  bsTooltip(id = 'min_peakwidth', title = 'Minimum peak width in seconds', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = 20, min = 0),
                  bsTooltip(id = 'max_peakwidth', title = 'Maximum peak width in seconds', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'snthresh', label = 'Signal to noise threshold', value = 10, min = 0),
                  bsTooltip(id = 'snthresh', title = 'Signal to noise ratio cut-off', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'prefilter', label = 'Prefilter', value = 3, min = 0),
                  bsTooltip(id = 'prefilter', title = 'Prefilter', placement = 'left', trigger = 'hover'),
                  numericInput(inputId = 'v_prefilter', label = 'Value of prefilter', value = 100, min = 0),
                  bsTooltip(id = 'v_prefilter', title = 'Value of prefilter', placement = 'left', trigger = 'hover'),
                  actionButton(inputId = 'peakdetectrun', label = 'Perform peak detection'),
                  actionButton(inputId = 'paramdetectrun', label = 'Detect parameters automatically')
                  ),
                div(id = 'align_param_box',
                    box(width = 2,
                        title = 'Input parameters for peak alignment',
                        selectInput(inputId = 'rtmethod', label = 'Method', choices = list('loess' = 'loess', 'obiwarp' = 'obiwarp'),selected = 'loess'),
                        # loess
                        div(id = 'loessparams',
                            numericInput(inputId = 'extra', label = 'Extra', value = 1, min = 0),
                            bsTooltip(id = 'extra', title = 'Extra', placement = 'left', trigger = 'hover'),
                            numericInput(inputId = 'span', label = 'Span', value = 0.25, min = 0),
                            bsTooltip(id = 'span', title = 'Span', placement = 'left', trigger = 'hover')),
                        # obiwarp
                        div(id = 'obiwarpparams',
                            numericInput(inputId = 'prof_step', label = 'profStep', value = 100, min = 0),
                            bsTooltip(id = 'prof_step', title = 'Prof step', placement = 'left', trigger = 'hover')),
                        # Grouping:
                        numericInput(inputId = 'bw', label = 'Bandwith', value = 10, min = 1),
                        bsTooltip(id = 'bw', title = 'Bandwith', placement = 'left', trigger = 'hover'),
                        numericInput(inputId = 'min_fraction', label = 'minFraction', value = 0.5, min = 0),
                        bsTooltip(id = 'min_fraction', title = 'Minimal fraction of samples a feature has to be present in', placement = 'left', trigger = 'hover'),
                        numericInput(inputId = 'min_samples', label = 'minSamples', value = 1, min = 1),
                        bsTooltip(id = 'min_samples', title = 'Minimum amount of samples a feature has to be present in', placement = 'left', trigger = 'hover'),
                        selectInput(inputId = 'fitgauss', label = 'Fitgauss', choices = list('False' = 'FALSE', 'True' = 'TRUE'),selected = 'FALSE'),
                        bsTooltip(id = 'fitgauss', title = 'Fitgauss', placement = 'left', trigger = 'hover'),
                        selectInput(inputId = 'verbose_columns', label = 'Verbose columns', choices = list('False' = 'FALSE', 'True' = 'TRUE'),selected = 'FALSE'),
                        bsTooltip(id = 'verbose_columns', title = 'Verbose columns', placement = 'left', trigger = 'hover'),
                        numericInput(inputId = 'integrate', label = 'Integrate', value = 1, min = 0, max = 1),
                        bsTooltip(id = 'integrate', title = 'Integrate', placement = 'left', trigger = 'hover'),
                        numericInput(inputId = 'max_features', label = 'maxFeatures', value = 100, min = 0),
                        bsTooltip(id = 'max_features', title = 'Maximum amount of features', placement = 'left', trigger = 'hover')
                    )
                ),
                box(width = 2,
                  title = 'Upload or save parameters',
                  downloadButton(outputId = 'save_params', label = 'Save parameters'),
                  fileInput(inputId = 'upload_params', label = 'Upload parameters:', multiple = F, accept = '.RData')
                ),
                box(
                  # width = 6,
                  title = 'Results', 
                  textOutput(outputId = 'peakamount'),
                  imageOutput(outputId = 'foundpeaks')
                )
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
    left_text = 'MetabOracle 0.2',
    right_text = 'Made by Micah Mallée'
  )
)

server <- function(input, output, session){
  rvalues <- reactiveValues()

  
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
  

  
  observe({
    rvalues$param_initial <- rvalues$parameters()
    print(rvalues$param_initial$RT_method)
    updateNumericInput(session = session, inputId = 'ppm', label = 'ppm', value = rvalues$param_initial$ppm, min = 0)
    updateNumericInput(session = session, inputId = 'noise', label = 'Noise', value = rvalues$param_initial$noise, min = 0)
    updateNumericInput(session = session, inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = rvalues$param_initial$min_peakwidth, min = 0)
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
  
  output$file <- renderTable(input$data_input)
  
  output$paths <- renderText(length(input$data_input$datapath))
  
  observeEvent(input$run, {
    if(input$inspect_trim == 1){
      output$inspect_plot <- renderPlot(PerformDataInspect(input$data_input$datapath))
      }
    else {
      rvalues$raw_data <- readMSData(files = input$data_input$datapath, mode = 'onDisk')
      output$inspect_plot <- renderPlot(plot(chromatogram(rvalues$raw_data)))
      }
    })
  
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
  
  observeEvent(input$peakdetectrun, {
    params <- SetPeakParam(platform = 'general', Peak_method = 'centWave', ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
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
      mSet <- PerformPeakAlignment(mSet, param = updateRawSpectraParam(param_optimized))
      mSet <- PerformPeakFiling(mSet, param = updateRawSpectraParam(param_optimized))
    }
    output$peakamount <- renderText(mSet[["msFeatureData"]][["chromPeakData"]]@nrows)
    chr <- chromatogram(mSet[['onDiskData']])
    xchr <- as(chr, 'XChromatograms')
    xchr[[1]]@chromPeaks <- mSet[["msFeatureData"]][["chromPeaks"]]
    xchr[[1]]@chromPeakData <- mSet[["msFeatureData"]][["chromPeakData"]]
    output$foundpeaks <- renderPlot(plot(xchr[[1]]))
  })
  
  observeEvent(input$paramdetectrun, {
    param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = 'loess', ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
    rvalues$optimized_params <- PerformParamsOptimization(raw_data = subset_raw_data, param = param_initial, ncore = 8)
    updateNumericInput(session = session, inputId = 'ppm', label = 'ppm', value = rvalues$optimized_params$best_parameters$ppm , min = 0)
    updateNumericInput(session = session, inputId = 'noise', label = 'Noise', value = rvalues$optimized_params$best_parameters$noise, min = 0)
    updateNumericInput(session = session, inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = rvalues$optimized_params$best_parameters$min_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = rvalues$optimized_params$best_parameters$max_peakwidth, min = 0)
    updateNumericInput(session = session, inputId = 'snthresh', label = 'Signal to noise threshold', value = rvalues$optimized_params$best_parameters$snthresh, min = 0)
    updateNumericInput(session = session, inputId = 'prefilter', label = 'Prefilter', value = rvalues$optimized_params$best_parameters$prefilter, min = 0)
    updateNumericInput(session = session, inputId = 'v_prefilter', label = 'Value of prefilter', value = rvalues$optimized_params$best_parameters$value_of_prefilter, min = 0)
  })
  
  output$save_params <- downloadHandler(
    filename = function() {
      paste('params_', Sys.Date(), '.RData', sep = '')
    },
    content = function(file) {
      param_initial <- SetPeakParam(platform = 'general', Peak_method = 'centWave', RT_method = input$rtmethod, ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                                            snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
      save(param_initial, file = file)
    }
  )
}

shinyApp(ui, server)

