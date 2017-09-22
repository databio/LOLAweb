options(shiny.maxRequestSize=100*1024^2)
source("themes.R")

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)

# load region data
dbPath = system.file("extdata", "hg19", package="LOLA")
regionDB = loadRegionDB(dbLocation=dbPath)

ui <- fluidPage(
  
  titlePanel("LOLA"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("userset", "Upload User Set(s)",
                multiple = TRUE,
                accept = c(".bed")),
      checkboxInput("checkbox", label = "Check Here to Upload Your Own Universe", value = FALSE),
      uiOutput("universe"),
      actionButton("run","runLOLA")
    ),
    mainPanel(
      plotOutput("plot"),
      DT::dataTableOutput("res")
    )
  )
)

server <- function(input, output) {
    
    output$universe <- renderUI({
      
      if(input$checkbox) {
        
        fileInput("useruniverse", "Upload Universe",
                  accept = c(".bed"))
        
      } else {
        
        selectInput("defaultuniverse", 
                    label = "Select Pre-Loaded Universe", 
                    choices = list.files("universes"))
        
      }
      
    })
    
    dat <- eventReactive(input$run, {
      
      withProgress(message = 'reticulating splines ... ', style = "old", value = 0, {
        
      userSets <- list()
      
      for (i in 1:length(input$userset[,1])) {
        
        userSet <- read.table(input$userset[[i, 'datapath']], header = F)
        colnames(userSet) <- c('chr','start','end','id','score','strand')
        userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
        
        userSets[[i]] <- userSet
      
      }
      
      userSets = GRangesList(userSets)
      
      if(input$checkbox) {
        
        userUniverse = read.table(file = input$useruniverse$datapath, header = F)
        
      } else {
        
        datapath <- paste0("universes/", input$defaultuniverse)
        
        userUniverse = read.table(file = datapath, header = F)
        
      }
      colnames(userUniverse) <- c('chr','start','end','id','score','strand')
      userUniverse <- with(userUniverse, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
      
      userSetsRedefined =	redefineUserSets(userSets, userUniverse)
      
      resRedefined = runLOLA(userSetsRedefined, userUniverse, regionDB, cores=2)
      
      return(resRedefined)
      
    })
  })
  
  output$plot <- renderPlot({
    
    ggplot(dat(), aes(logOddsRatio)) +
      geom_histogram() +
      theme_ns()

  })
  
  output$res <- DT::renderDataTable({
  
    dat()
    
  })
}

shinyApp(ui = ui, server = server)
