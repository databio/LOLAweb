options(shiny.maxRequestSize=100*1024^2)
source("themes.R")
source("misc.R")
source("disabler.R")

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)
library(DT)

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  titlePanel("LOLA"),
  
      fluidRow(
        column(8,
               fileInput("userset", "Upload User Set(s)",
                         multiple = TRUE,
                         accept = c(".bed")),
               tags$div(
                h3("Select a universe",
                tags$a(href = "http://code.databio.org/LOLA/articles/choosingUniverse.html", 
                       icon("question-circle-o"), 
                       target = "blank"))
               ),
               uiOutput("universe"),
               checkboxInput("checkbox", 
                             label = "Check Here to Upload Your Own Universe",
                             value = FALSE),
               tags$div(
               h3("Select Region Database",
               tags$a(href = "http://databio.org/regiondb", 
                      icon("question-circle-o"), 
                      target = "blank"))
               ),
               # radioButtons("loladb", 
               #              "", 
               #              choices = c("Core", "Extended")), 
               # must be a better way to do this
               HTML(disabledbutton),
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
                            class = "runLOLA"),
               textOutput("messages")
        ),
      class = "headerBox"),
      fluidRow(
        
        column(2,
                 uiOutput("slider_logOdds"),
                 uiOutput("slider_support"),
                 uiOutput("slider_pvalue"),
                 uiOutput("select_collection"),
                 uiOutput("select_sort"),
                 uiOutput("select_userset")
                 
          ),

        column(5,
               conditionalPanel(condition = "output.res",
                                h3("Log Odds"),
                                downloadButton("logodds_plot_dl",
                                               label = "Download Plot")),
               plotOutput("logodds_plot"),
               conditionalPanel(condition = "output.res",
                                h3("Support"),
                                downloadButton("support_plot_dl",
                                               label = "Download Plot")),
               plotOutput("support_plot")
        ),
        column(5,
               conditionalPanel(condition = "output.res",
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
    
    
  observeEvent(input$run, {
    withCallingHandlers({
      shinyjs::html(id = "messages", html = "")
      raw_dat()
    },
    message = function(m) {
      shinyjs::html(id = "messages", html = m$message, add = FALSE)
    },
    warning = function(m) {
      shinyjs::html(id = "messages", html = m$message, add = FALSE)
    })
  })
  
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
    
    raw_dat <- eventReactive(input$run, {
      
      withProgress(message = 'calculating region set enrichments...', style = "old", value = 0, {
        
      userSets <- list()
      
      for (i in 1:length(input$userset[,1])) {
        
        userSet <- read.table(input$userset[[i, 'datapath']], header = F)
        colnames(userSet) <- c('chr','start','end','id','score','strand')
        userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
        
        userSets[[i]] <- userSet
      
      }
      
      userSets = GRangesList(userSets)
      
      names(userSets) = input$userset[,"name"]
      
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
      
      # cores = parallel::detectCores()
      
      resRedefined = runLOLA(userSetsRedefined, 
                             userUniverse, 
                             regionDB,
                             cores=1)
      
      # need to make sure user set is discrete even if coded as number
      resRedefined$userSet = as.character(resRedefined$userSet)
      
      # seeing some missing pvalues need to make sure these are not present
      # resRedefined = subset(resRedefined, 
      #                       logOddsRatio != "",
      #                       pValueLog != "",
      #                       support != ""
      #                       )
      
      return(resRedefined)
      
    })
  })
    
  dat <- reactive({
    
    dat <- subset(raw_dat(), 
                  logOddsRatio >= input$slider_logOdds_i & 
                    support >= input$slider_support_i &
                    pValueLog >= input$slider_pvalue_i)
  
    if (input$select_collection_i != "All Collections") {
      
       dat <- subset(dat, collection == input$select_collection_i)
       
    }
    
    if (input$select_userset_i != "All User Sets") {

      dat <- subset(dat, userSet == input$select_userset_i)

    }
    
    dat$id <- paste(dat$description, dat$dbSet, sep = "_")
      
    return(dat)

  })
    
  setchoices <- function() {
    
    req(input$run)
    
    if(length(unique(raw_dat()$userSet)) == 1) {
      
      unique(raw_dat()$userSet)
      
    } else {
      
      c("All User Sets", unique(raw_dat()$userSet))
      
    }
  }
  
  output$select_collection <- renderUI({
    
    req(input$run)
    
    list(
      HTML("<hr>"),
      selectInput("select_collection_i", 
                "Select Collection", 
                choices = c("All Collections", unique(raw_dat()$collection)),
                selected = "All Collections"))
    
  })  
  
  output$select_sort <- renderUI({
    
    req(input$run)
    
    selectInput("select_sort_i", 
                "Select Sort Column", 
                choices = names(raw_dat()),
                selected = "meanRank")
    
  })  
  
  output$select_userset <- renderUI({
    
    req(input$run)
    
    selectInput("select_userset_i", 
                "Select User Set", 
                choices = setchoices())
    
  })  
  
  # barplots
    
  # log odds
  
  output$slider_logOdds <- renderUI({
    
    req(input$run)
    
    sliderInput("slider_logOdds_i", 
                "Specify Log Odds Cutoff", 
                min = round(min(raw_dat()$logOddsRatio), 3), 
                max = round(max(raw_dat()$logOddsRatio), 3),
                value = round_top(raw_dat()$logOddsRatio, 30))
    
    })  
  
  # set up function
  logodds_plot_input <- function() {
    
    ggplot(dat(), aes(reorder(description, eval(parse(text = input$select_sort_i))), logOddsRatio, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("Log Odds Ratio") +
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
  
  # support
  
  # slider
  output$slider_support <- renderUI({
    
    req(input$run)
    
    sliderInput("slider_support_i", 
                "Specify Support Cutoff", 
                min = round(min(raw_dat()$support), 3), 
                max = round(max(raw_dat()$support), 3),
                value = round(min(raw_dat()$support), 3))
    
  })  
  
  # set up function
  support_plot_input <- function() {
    
    ggplot(dat(), aes(reorder(description, eval(parse(text = input$select_sort_i))), support, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("Support") +
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
  
  # pvalue
  
  # slider input
  
  output$slider_pvalue <- renderUI({
    
    req(input$run)
    
    sliderInput("slider_pvalue_i", 
                "Specify P Value Cutoff", 
                min = round(min(raw_dat()$pValueLog), 3), 
                max = round(max(raw_dat()$pValueLog), 3),
                value = round_top(raw_dat()$pValueLog, 30))
    
  })  
  
  # set up function
  
  pvalue_plot_input <- function() {
    
    ggplot(dat(), aes(reorder(description, eval(parse(text = input$select_sort_i))), pValueLog, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("P Value (log scale)") +
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
    
    dat <- dat()[order(eval(parse(text = input$select_sort_i)), decreasing = TRUE)]
    
    dat$dbSet <- ifelse(dat$collection == "sheffield_dnase",
                        paste0("<a href = 'http://db.databio.org/clusterDetail.php?clusterID=",
                               tools::file_path_sans_ext(dat$filename),
                               "' target = '_blank'>",
                                dat$dbSet,
                                "</a>"),
                        dat$dbSet)
    
    dat %>%
      datatable(rownames = FALSE, 
                escape = FALSE,
                extensions = c("Responsive", "Buttons"),
                options = list(dom = "Bfrtip",
                               buttons = "csv")) %>%
      formatRound(columns=c('logOddsRatio', 'pValueLog'),
                  digits = 4)
  })
  
  # observe({
  #   
  #   if(nrow(dat()) == 0) {
  #     showNotification("No results found. Please try adjusting the inputs.", 
  #                      type = "error", 
  #                      duration = 10)
  #   } else {
  #     
  #     invisible()
  #     
  #   }
  #   
  # })
  
}

shinyApp(ui = ui, server = server)
