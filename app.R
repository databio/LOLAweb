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
  
  titlePanel(title = HTML("<img src='LOLA-logo.png' alt='LOLA logo' width='200'>"),
             windowTitle = "LOLA"),
  
      fluidRow(
        column(4,
               tags$div(
                 h3("#1 Select User Set",
                    tags$a(href = "http://code.databio.org/LOLA/articles/gettingStarted.html", 
                           icon("question-circle-o"), 
                           target = "blank"))
               ),
               uiOutput("usersets"),
               checkboxInput("checkbox_userset", 
                             label = "Check Here to Upload Your Own User Set(s)",
                             value = TRUE)
        ),
        column(4,
               tags$div(
                 h3("#2 Select Universe",
                    tags$a(href = "http://code.databio.org/LOLA/articles/choosingUniverse.html", 
                           icon("question-circle-o"), 
                           target = "blank"))
               ),
               uiOutput("universe"),
               checkboxInput("checkbox", 
                             label = "Check Here to Upload Your Own Universe",
                             value = FALSE),
               actionButton("run",
                            "RUN LOLA", 
                            class = "runLOLA"),
               HTML("<br>"),
               HTML("<br>"),
               column(1,
                      htmlOutput("gear")
                      ),
               column(3,
                      htmlOutput("messages")
                      )
               ),
        column(4, 
               tags$div(
                 h3("#3 Select Region Database",
                    tags$a(href = "http://databio.org/regiondb", 
                           icon("question-circle-o"), 
                           target = "blank"))
               ),
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
      class = "headerBox"),
      fluidRow(
        
        column(2,
                 uiOutput("slider_rank"),
                 uiOutput("slider_oddsratio"),
                 uiOutput("slider_support"),
                 uiOutput("slider_pvalue"),
                 uiOutput("select_collection"),
                 uiOutput("select_sort"),
                 uiOutput("select_userset")
                 
          ),

        column(5,
               conditionalPanel(condition = "output.res",
                                h3("Odds Ratio"),
                                downloadButton("oddsratio_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("oddsratio_plot"),
               conditionalPanel(condition = "output.res",
                                h3("Support"),
                                downloadButton("support_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("support_plot")
        ),
        column(5,
               conditionalPanel(condition = "output.res",
                                h3("P Value"),
                                downloadButton("pvalue_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
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
        shinyjs::html(id = "gear", html = "<i class='fa fa-2x fa-spin fa-cog'></i>", add = FALSE)
        raw_dat()
      },
      message = function(m) {
        shinyjs::html(id = "messages", html = m$message, add = FALSE)
        
      })
    
      shinyjs::hide("messages")
      shinyjs::hide("gear")
      
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
    
    output$usersets <- renderUI({
      
      if(input$checkbox_userset) {
        
        fileInput("userset", "Upload User Set(s)",
                  multiple = TRUE,
                  accept = c(".bed"))
        
      } else {
        
        selectInput("defaultuserset", 
                    label = "Select Pre-Loaded User Set", 
                    choices = list.files("userSets"))
        
      }
      
    })
    
    raw_dat <- eventReactive(input$run, {
      
      message("Calculating region set enrichments ...")
      userSets <- list()
      
      if(input$checkbox_userset) {
        
        for (i in 1:length(input$userset[,1])) {
          
          userSet <- read.table(input$userset[[i, 'datapath']], header = F)
          colnames(userSet) <- c('chr','start','end','id','score','strand')
          userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
          
          userSets[[i]] <- userSet
          
        }
        
        userSets = GRangesList(userSets)
        
        names(userSets) = input$userset[,"name"]

      } else {
        
        datapath <- paste0("userSets/", input$defaultuserset)
        
        userSet = read.table(file = datapath, header = F)
        colnames(userSet) <- c('chr','start','end','id','score','strand')
        userSet <- with(userSet, GRanges(chr, IRanges(start+1, end), strand, score, id=id))
        
        userSets[[1]] <- userSet
        
        userSets = GRangesList(userSets)
        
        names(userSets) = input$defaultuserset
        
      }
      
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
      #                       oddsRatio != "" &
      #                       pValueLog != "" &
      #                       support != ""
      #                       )
      
      # resRedefined = subset(resRedefined,
      #                       oddsRatio > 0 &
      #                       pValueLog > 0 &
      #                       support > 0
      #                       )
      
      # garbage collection
      # rm(regionDB, userSetsRedefined, userUniverse)
      
      return(resRedefined)
      

  })
    
  dat <- reactive({
    
    dat <- subset(raw_dat(), maxRnk <= input$slider_rank_i)
    
    dat <- subset(dat,
                  oddsRatio >= input$slider_oddsratio_i &
                    support >= input$slider_support_i &
                    pValueLog >= input$slider_pvalue_i)
    
    # dat <- subset(dat, maxRnk >= input$slider_rank_i)
  
    if (input$select_collection_i != "All Collections") {
      
       dat <- subset(dat, collection == input$select_collection_i)
       
    }
    
    if (input$select_userset_i != "All User Sets") {

      dat <- subset(dat, userSet == input$select_userset_i)

    }
    
    dat$id <- paste(dat$description, dat$dbSet, sep = "_")
    
    dat$axis_label <- strtrim(dat$description, 50)
    
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
  
  output$slider_rank <- renderUI({
    
    req(input$run)
    
    sliderInput("slider_rank_i", 
                "Specify Max Rank Cutoff", 
                min = min(raw_dat()$maxRnk),
                max = max(raw_dat()$maxRnk),
                value = quantile(raw_dat()$maxRnk, probs = 20/nrow(raw_dat())))
    
  })  
  
  output$select_sort <- renderUI({
    
    req(input$run)
    
    selectInput("select_sort_i", 
                "Select Sort Column", 
                choices = names(raw_dat()),
                selected = "meanRnk")
    
  })  
  
  output$select_userset <- renderUI({
    
    req(input$run)
    
    selectInput("select_userset_i", 
                "Select User Set", 
                choices = setchoices())
    
  })  
  
  # barplots
    
  # odds ratio
  
  output$slider_oddsratio <- renderUI({
    
    req(input$run)
    
    sliderInput("slider_oddsratio_i",
                "Specify Odds Ratio Cutoff",
                min = round(min(raw_dat()$oddsRatio), 3),
                max = round(max(raw_dat()$oddsRatio), 3),
                value = round(min(raw_dat()$oddsRatio), 3))
    })
  
  # set up function
  oddsratio_plot_input <- function() {
    
    ggplot(dat(), aes(reorder(axis_label, eval(parse(text = input$select_sort_i))), oddsRatio, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("Odds Ratio") +
      coord_flip() +
      theme_ns()
    
  }
  
  # call plot
  output$oddsratio_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    oddsratio_plot_input()

  })
  
  # download handler
  output$oddsratio_plot_dl <- downloadHandler(
    filename = function() { paste(gsub(".bed", "", input$userset),
                                  "_oddsratio", 
                                  ".png", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = oddsratio_plot_input(), device = "png")
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
    
    ggplot(dat(), aes(reorder(axis_label, eval(parse(text = input$select_sort_i))), support, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("Support") +
      coord_flip() +
      theme_ns()
    
  }
  
  output$support_plot <- renderPlot({
    
    req(input$select_sort_i)
    
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
                value = round(min(raw_dat()$pValueLog), 3))
    
    
  })  
  
  # set up function
  
  pvalue_plot_input <- function() {
    
    ggplot(dat(), aes(reorder(axis_label, eval(parse(text = input$select_sort_i))), pValueLog, fill = userSet, group = id)) +
      geom_bar(stat = "identity", position = "dodge") +
      xlab("Description") +
      ylab("P Value (log scale)") +
      coord_flip() +
      theme_ns()
    
  }
  output$pvalue_plot <- renderPlot({
    
    req(input$select_sort_i)
    
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
    
    req(input$select_sort_i)
    
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
                               buttons = list(list(
                                 extend = "csv",
                                 text = '<i class="fa fa-download"></i> Download CSV')),
                               paging = FALSE)) %>%
      formatRound(columns=c('oddsRatio', 'pValueLog'),
                  digits = 4)
  })
  
}

shinyApp(ui = ui, server = server)
