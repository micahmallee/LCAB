library(shiny)
library(MetaboAnalystR)
library(shinyFiles)

options(shiny.maxRequestSize=100*1024^2)

# Define UI
ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      helpText("Welcome to MetabOracle! Please upload your mzxml data below."),
      fileInput(inputId = "raw_data", label = "Choose mzXML file", multiple = TRUE, accept = '.mzxml',),
      sliderInput(inputId = 'rt.idx', label = 'rt.idx', min = 0, max = 1, step = 0.1, value = 0.6),
      radioButtons(inputId = 'inspect_trim', label = 'Inspect or trim data.',
                         choices = list('Inspect' = 1, 
                                        'Trim' = 2), 
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
  output$file <- renderTable(input$raw_data
  )
  
  output$type <- renderText(input$raw_data$datapath)
  
  observeEvent(input$run, {
    if(input$inspect_trim == 1){
      output$inspect_plot <- renderPlot(PerformDataInspect(input$raw_data$datapath))
    }
    else {
      trimmed_data <- PerformDataTrimming(datapath = input$raw_data$datapath, rt.idx = input$rt.idx)
      output$inspect_plot <- renderPlot(plot(chromatogram(trimmed_data)))
    }
  })
}

shinyApp(ui, server)

