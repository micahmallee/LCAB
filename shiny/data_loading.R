library(shiny)
library(MetaboAnalystR)
library(shinyFiles)
library(bs4Dash)
library(shinydashboard)
library(shinydashboardPlus)

options(shiny.maxRequestSize=100*1024^2)

# Define UI
ui <- dashboardPagePlus(
  skin = 'midnight',
  dashboardHeaderPlus(
    title = tagList(
      tags$span(
        class = "logo-lg", "MetabOralce"
      )
    ),
    left_menu = tagList(uiOutput("admincontrol")),
    enable_rightsidebar = TRUE,
    rightSidebarIcon = "file-upload",
    userOutput("usertop")
  ),
  dashboardSidebar(
    collapsed = TRUE,
    sidebarMenu(
      menuItem('Dashboard', tabName = 'dashboard', icon = icon('dashboard')),
      menuItem('Peak Detection', tabName = 'peakpicking', icon = icon('search'))
    )
  ),
  
  rightsidebar = rightSidebar(collapsed = FALSE,
    uiOutput("sidebarright"),
    title = 'Load data',
    helpText("Welcome to MetabOracle! Please upload your mzXML/mzML/netCDF data below."),
    fileInput(inputId = "data_input", label = "Choose data file(s)"),
    sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
    radioButtons(inputId = 'inspect_trim', label = 'Inspect or Upload data.',
                 choices = list('Inspect' = 1,
                                'Upload' = 2),
                 selected = 1, inline = T),
    actionButton(inputId = 'run', label = 'Run')
  ),
  
  dashboardBody(
    tabItems(
      tabItem('dashboard', 
              fluidRow(
                box(
                  title = 'Output',
                  tableOutput(outputId = 'file'),
                  textOutput(outputId = 'type'),
                  textOutput(outputId = 'test')
                  ),
                box(
                  imageOutput(outputId = 'inspect_plot')
                  ))),
      tabItem('peakpicking',
              fluidRow(
                box(
                  title = 'Input parameters for CentWave algorithm',
                  numericInput(inputId = 'ppm', label = 'ppm', value = 5, min = 0),
                  numericInput(inputId = 'noise', label = 'Noise', value = 1000, min = 0),
                  numericInput(inputId = 'min_peakwidth', label = 'Minimal peakwidth', value = 5, min = 0),
                  numericInput(inputId = 'max_peakwidth', label = 'Maximum peakwidth', value = 20, min = 0),
                  numericInput(inputId = 'snthresh', label = 'Signal to noise threshold', value = 10, min = 0),
                  numericInput(inputId = 'prefilter', label = 'Prefilter', value = 3, min = 0),
                  numericInput(inputId = 'v_prefilter', label = 'Value of prefilter', value = 100, min = 0),
                  actionButton(inputId = 'peakdetectrun', label = 'Perform peak detection')
                  ),
                box(
                  imageOutput(outputId = 'inspect_plot'),
                  imageOutput(outputId = 'foundpeaks')
                ))
              )
    )
  )
)

server <- function(input, output){
  output$file <- renderTable(input$data_input
  )
  
  output$type <- renderText(input$data_input$datapath)
  
  observeEvent(input$run, {
    if(input$inspect_trim == 1){
      output$inspect_plot <- renderPlot(PerformDataInspect(input$data_input$datapath))
    }
    else {
      # oke <- input$data_input$datapath[[1]][-5]
      # oke <- paste(oke, collapse = '/')
      raw_data <- readMSData(files = input$data_input$datapath, mode = 'onDisk')
      output$inspect_plot <- renderPlot(plot(chromatogram(raw_data)))
    }
  })
  
  observeEvent(input$peakdetectrun, {
    params <- SetPeakParam(platform = 'general', Peak_method = 'centWave', ppm = input$ppm, noise = input$noise, min_peakwidth = input$min_peakwidth, max_peakwidth = input$max_peakwidth,
                           snthresh = input$snthresh, prefilter = input$prefilter, value_of_prefilter = input$v_prefilter)
    mSet <- PerformPeakPicking(raw_data, updateRawSpectraParam(params))
    chr <- chromatogram(mSet[['onDiskData']])
    xchr <- as(chr, 'XChromatograms')
    xchr[[1]]@chromPeaks <- mSet[["msFeatureData"]][["chromPeaks"]]
    xchr[[1]]@chromPeakData <- mSet[["msFeatureData"]][["chromPeakData"]]
    output$foundpeaks <- renderPlot(plot(xchr))
  })
}

shinyApp(ui, server)

