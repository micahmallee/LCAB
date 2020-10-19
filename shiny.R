library(shiny)
library(MetaboAnalystR)
library(xcms)


# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Upload mzXML"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose mzXML File",
                multiple = FALSE,
                accept = c(".mzXML")),
      
      # Horizontal line ----
      tags$hr(),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents")
      
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderText({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        data <- PerformDataInspect(datapath = input)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(input$disp == "head") {
      return(head(data))
    }
    else {
      return(data)
    }
    
  })
  
}


shinyApp(ui, server)
