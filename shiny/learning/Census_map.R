library(shiny)
library(maps)
library(mapproj)
source('shiny/learning/helpers.R')
counties <- readRDS("shiny/learning/counties.rds")

# Define UI
ui <- fluidPage(
  titlePanel('Prototype files inladen'),
  
  navbarPage("MetabOracle",
             navbarMenu('Preprocessing',
                        tabPanel("Load data"),
                        tabPanel('Inspect data')
                        )
             ),
  
  sidebarLayout(
    position = 'left',
    sidebarPanel(
      helpText('Welcome to MetabOracle, here you can analyse all your metabolomic data!'),
      selectInput("var",
                  label = "Choose variable to display",
                  choices = c("Percent White",
                              "Percent Black",
                              "Percent Hispanic",
                              "Percent Asian"),
                  selected = "Percent White"),
      sliderInput("range",
                  label = "Range of interest:",
                  min = 0, max = 100, value = c(0, 100))
      ),
    
    mainPanel(
      p('Made mainly for the Rhinoceromics project.'),
      textOutput("selected_var"),
      textOutput("min_max"),
      plotOutput("map"))
  )
  
)

# Define server logic
server <- function(input, output) {
  output$map <- renderPlot({
    args <- switch(input$var,
                   "Percent White" = list(counties$white, "darkgreen", "% White"),
                   "Percent Black" = list(counties$black, "black", "% Black"),
                   "Percent Hispanic" = list(counties$hispanic, "darkorange", "% Hispanic"),
                   "Percent Asian" = list(counties$asian, "darkviolet", "% Asian"))
    
    args$min <- input$range[1]
    args$max <- input$range[2]
    
    do.call(percent_map, args)
  })
}


# Run the app
shinyApp(ui, server)
