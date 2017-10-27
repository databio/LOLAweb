options(shiny.maxRequestSize=100*1024^2)
source("themes.R")

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)
library(DT)

ui <- fluidPage(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  titlePanel("LOLA"),
  
      fluidRow(
        column(8,
               fileInput("userset", "Upload User Set(s)",
                         multiple = TRUE,
                         accept = c(".bed")),
               checkboxInput("checkbox", 
                             label = "Check Here to Upload Your Own Universe",
                             value = FALSE),
               uiOutput("universe"),
               radioButtons("loladb", 
"Specify Core or Extended LOLA Databases", choices = c("Core", "Extended")), 
               conditionalPanel(condition = "input.loladb == 'Core'",
                                selectInput("refgenome_core", 
                                            "Reference Genome", 
                                            choices = c("hg19","hg38", "mm10"))),
               conditionalPanel(condition = "input.loladb == 'Extended'",
                                selectInput("refgenome_ext", 
                                            "Reference Genome", 
                                            choices = c("hg19","hg38")))
        ),
        column(4, 
               actionButton("run",
                            "RUN LOLA", 
                            class = "runLOLA")
        ),
      class = "headerBox"),
      fluidRow(
        column(4,
               conditionalPanel(condition = "output.logodds_plot",
                                h3("Log Odds"),
                                downloadButton("logodds_plot_dl",
                                               label = "Download Plot")),
               plotOutput("logodds_plot")
        ),
        column(4, 
               conditionalPanel(condition = "output.support_plot",
                                h3("Support"),
                                downloadButton("support_plot_dl",
                                               label = "Download Plot")),
               plotOutput("support_plot")
        ),
        column(4,
               conditionalPanel(condition = "output.pvalue_plot",
                                h3("P Value"),
                                downloadButton("pvalue_plot_dl", 
                                               label = "Download Plot")),
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
      
      # load region data for each reference genome
      dbPath = paste("reference", 
                      input$loladb, 
                      ifelse(input$loladb == "Core", 
                             input$refgenome_core,
                             input$refgenome_ext),
                     sep = "/"
                      )
                      
      regionDB = loadRegionDB(dbLocation=dbPath)
      
      resRedefined = runLOLA(userSetsRedefined, 
                             userUniverse, 
                             regionDB,
                             cores=2)
      
      return(resRedefined)
      
    })
  })
  
  # barplots
    
  # log odds
    
  # set up function
  logodds_plot_input <- function() {
    
    ggplot(dat(), aes(description, logOddsRatio)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()
    
  }
  
  # call plot
  output$logodds_plot <- renderPlot({
    
    logodds_plot_input()

  })
  
  # download handler
  output$logodds_plot_dl <- downloadHandler(
    filename = function() { paste(gsub(".bed", "", input$userset),
                                  "_logodds", 
                                  ".png", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = logodds_plot_input(), device = "png")
    }
  )
  
  # set up function
  support_plot_input <- function() {
    
    ggplot(dat(), aes(description, support)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()
    
  }
  
  output$support_plot <- renderPlot({
    
    support_plot_input()
    
  })
  
  # download handler
  output$support_plot_dl <- downloadHandler(
    filename = function() { paste(gsub(".bed", "", input$userset),
                                  "_support", 
                                  ".png", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = support_plot_input(), device = "png")
    }
  )
  
  # set up function
  
  pvalue_plot_input <- function() {
    
    ggplot(dat(), aes(description, pValueLog)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_ns()
    
  }
  output$pvalue_plot <- renderPlot({
    
    pvalue_plot_input()
    
  })
  
  output$pvalue_plot_dl <- downloadHandler(
    filename = function() { paste(gsub(".bed", "", input$userset),
                                  "_pvalue", 
                                  ".png", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = pvalue_plot_input(), device = "png")
    }
  )
  
  output$res <- DT::renderDataTable({
    
    # data.table of results ordered by ranks
    dat()[order(meanRnk, maxRnk)] %>%
      datatable(rownames = FALSE) %>%
      formatRound(columns=c('logOddsRatio', 'pValueLog'),
                  digits = 4)
    
  })
}

shinyApp(ui = ui, server = server)
