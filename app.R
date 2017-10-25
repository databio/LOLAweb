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
  
      fluidRow(
        column(8,
                   fileInput("userset", "Upload User Set(s)",
                             multiple = TRUE,
                             accept = c(".bed")),
                   checkboxInput("checkbox", label = "Check Here to Upload Your Own Universe", value = FALSE),
                   uiOutput("universe")
        ),
        column(4, 
            # actionButton("run","runLOLA", 
            #              style='padding:40px; font-size:200%')
            actionButton("run","runLOLA")
              
        )
      ),
      fluidRow(
        column(4,
               plotOutput("logodds_plot")
        ),
        column(4, 
               plotOutput("support_plot")
        ),
        column(4,
               plotOutput("pvalue_plot")
        )
      ),
      fluidRow(
        column(DT::dataTableOutput("res"), width = 12) 
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
      
      withProgress(message = 'calculating region set enrichments...', style = "old", value = 0, {
        
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
  
  output$logodds_plot <- renderPlot({
    
    ggplot(dat(), aes(description, logOddsRatio)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()

  })
  
  output$support_plot <- renderPlot({
    
    ggplot(dat(), aes(description, support)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()
    
  })
  
  output$pvalue_plot <- renderPlot({
    
    ggplot(dat(), aes(description, pValueLog)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()
    
  })
  
  output$res <- DT::renderDataTable({
    
    # data.table of results ordered by ranks
    dat()[order(meanRnk, maxRnk)]
    
  }, options = list(scrollX = TRUE))
}

shinyApp(ui = ui, server = server)
