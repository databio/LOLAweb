library(shiny)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("PROD"),
  
  mainPanel(
    sliderInput("duration", "Duration (minutes)", min = 1, max = 60, value = 15),
    actionButton("start", "Start"),
    actionButton("stop", "Stop"),
    actionButton("reset", "Reset"),
    textOutput("countdown")
  )
)

server <- function(input, output,session) {
  
  running <- reactiveVal(TRUE)
  # # observers for actionbuttons
  # observeEvent(input$start, {active(TRUE)})
  # observeEvent(input$stop, {active(FALSE)})
  # observeEvent(input$reset, {timer(input$duration*60)})
  # 
  # # Initialize the timer, 10 seconds, not active.
  # 
  # timer <- reactiveVal(1000)
  # active <- reactiveVal(FALSE)
  # 
  # # Output the time left.
  # output$countdown <- renderText({
  #   
  #   req(input$start)
  #   
  #   dur <- timer()
  #   
  #   mins <- format(floor(dur/60), digits = 2)
  #   secs <- format(dur %% 60, digits = 2)
  #   
  #   paste0(mins, ":",secs)
  #   
  # })
  # 
  # # observer that invalidates every second. If timer is active, decrease by one.
  # observe({
  #   invalidateLater(1000, session)
  #   isolate({
  #     if(active())
  #     {
  #       timer(timer() - 1)
  #       if(timer() < 1)
  #       {
  #         active(FALSE)
  #         showModal(modalDialog(
  #           title = "Coundown complete",
  #           "Time's up!"
  #         ))
  #       }
  #     }
  #   })
  # })
  
  clock_time <- eventReactive(input$start, {

    dur <- input$duration*60

    while(dur > 0 & isolate(running())) {

      Sys.sleep(0.99)

      dur <- dur - 1

      mins <- format(floor(dur/60), digits = 2)
      secs <- format(dur %% 60, digits = 2)

      message(paste0(mins, ":",secs))

    }

    showModal(modalDialog("Time's up!"))


  })


  observeEvent(input$start, {

    withCallingHandlers({
      shinyjs::html(id = "countdown", html = "")
      # shinyjs::html(id = "gear", html = "<i class='fa fa-4x fa-spin fa-cog'></i>", add = FALSE)
      clock_time()
    },
    message = function(m) {
      shinyjs::html(id = "countdown", html = m$message, add = FALSE)

    })

  })
  
  observeEvent(input$stop, {
    
    running <- FALSE
    
  })

  
}

# Run the application 
shinyApp(ui = ui, server = server)

