library(shiny)
library(MetaboAnalystR)

# Define UI
ui <- fluidPage(
  
)

server <- function(input, output){
  output$inspection <- renderPlot(
    args <- switch(input$var,
                   )
  )
  
  do.call(MetaboAnalystR::PerformDataInspect(args))
}

shinyApp(ui, server)