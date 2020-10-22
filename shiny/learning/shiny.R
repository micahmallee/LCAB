library(shiny)

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
      textOutput("min_max"))
  )
  
)

# Define server logic
server <- function(input, output) {
  output$selected_var <- renderText({
    paste("You have selected:", input$var)
  })
  
  output$min_max <- renderText({
    paste("You have selected a range that goes from", input$range[1], "to", input$range[2])
  })
}


# Run the app
shinyApp(ui, server)
