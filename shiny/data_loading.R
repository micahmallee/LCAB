library(shiny)
library(MetaboAnalystR)
library(shinyFiles)

options(shiny.maxRequestSize=100*1024^2)

# Define UI
ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      helpText("Welcome to MetabOracle! Please upload your mzxml data below."),
      fileInput(inputId = "data_input", label = "Choose mzXML file", multiple = TRUE, accept = '.mzxml',),
      sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
      radioButtons(inputId = 'inspect_trim', label = 'Inspect or Upload data.',
                         choices = list('Inspect' = 1, 
                                        'Upload' = 2), 
                         selected = 1, inline = T),
      actionButton(inputId = 'run', label = 'Run')
      ),
    
    mainPanel(
      tableOutput(outputId = 'file'),
      textOutput(outputId = 'type'),
      imageOutput(outputId = 'inspect_plot'),
      textOutput(outputId = 'test')
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
      raw_data <- ImportRawMSData(foldername = input$data_input$datapath, plotSettings = SetPlotParam(Plot = F), mode = 'inMemory')
      output$inspect_plot <- renderPlot(plot(chromatogram(raw_data)))
    }
  })
}

shinyApp(ui, server)

