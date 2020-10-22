library(shiny)
library(MetaboAnalystR)

# Define UI
ui <- fluidPage(
  fileInput(inputId = "raw_data", label = "Choose mzXML file", multiple = TRUE, accept = '.mzxml',),
  mainPanel(
    textOutput('path')
  )
)

server <- function(input, output){
  output$path <- renderText(expr = input$raw_data
  )
}

shinyApp(ui, server)