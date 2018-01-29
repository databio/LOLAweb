options(shiny.maxRequestSize=100*1024^2)

source("themes.R")
source("misc.R")
source("disabler.R")

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)
library(DT)
library(simpleCache)
library(sodium)

setCacheDir("cache")

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
               shinyjs::hidden(
                 selectInput("defaultuserset", 
                                           label = "Select Pre-Loaded User Set",
                                           choices = list.files("userSets"))
                 ),
               fileInput("userset", "Upload User Set(s)",
                         multiple = TRUE,
                         accept = c(".bed")),
               actionButton("button_userset_example",
                            "Load example data"),
               shinyjs::hidden(
                 actionButton("button_userset_upload",
                            "Upload data")
               )
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
                            class = "runLOLA")
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
               htmlOutput("gear"),
               uiOutput("slider_rank"),
               uiOutput("slider_oddsratio"),
               uiOutput("slider_support"),
               uiOutput("slider_pvalue"),
               uiOutput("select_collection"),
               uiOutput("select_sort"),
               uiOutput("select_userset")
          ),
        column(10,
               htmlOutput("messages"),
               tags$h4(htmlOutput("link"))),
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

server <- function(input, output, session) {
  
  exampleuserset <- reactiveValues(toggle = TRUE)
  
  # get url parameter for retrieval query key
  query <- reactive({
    
    parseQueryString(session$clientData$url_search)
    
    })
  
  observeEvent(input$button_userset_upload, {
    
    list(shinyjs::show("button_userset_example"),
         shinyjs::hide("button_userset_upload"),
         shinyjs::hide("defaultuserset"),
         shinyjs::show("userset"),
         exampleuserset$toggle <- FALSE)
    
  })
  
  observeEvent(input$button_userset_example, {
    
    list(shinyjs::show("button_userset_upload"),
         shinyjs::hide("button_userset_example"),
         shinyjs::show("defaultuserset"),
         shinyjs::hide("userset"),
         exampleuserset$toggle <- TRUE)
    
  })
  
  observeEvent(input$run, {
      
    withCallingHandlers({
      shinyjs::html(id = "messages", html = "")
      shinyjs::html(id = "gear", html = "<i class='fa fa-4x fa-spin fa-cog'></i>", add = FALSE)
      raw_dat_nocache()
    },
    message = function(m) {
      shinyjs::html(id = "messages", html = m$message, add = FALSE)
        
    })
    
    shinyjs::hide("messages")
    shinyjs::hide("gear")
      
  })
    
  output$link <- renderText({
    
    req(input$select_sort_i)
    
    baseurl <- session$clientData$url_hostname
    port <- session$clientData$url_port
    pathname <- session$clientData$url_pathname

 
    link <- paste0(baseurl, ":", port, pathname, "?key=", keyphrase())
    
    paste0("<a href = 'http://", link, "' target = 'blank'>", link, "</a>")
    
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
    
  keyphrase <- eventReactive(input$run, {
      
    paste0(sample(c(LETTERS,1:9), 15), collapse = "")
      
    })
    
    raw_dat_nocache <- eventReactive(input$run, {
    
      message("Calculating region set enrichments ...")
      userSets <- list()

      # if(!input$switch_userset) {
      if(!exampleuserset$toggle) {


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

      # caching
      # need to call keyphrase from reactive above because it is used to construct link
      key <- hash(charToRaw(keyphrase()))
      msg <- serialize(resRedefined, connection = NULL)

      cipher <- data_encrypt(msg, key)

      simpleCache(keyphrase(), cipher)

      return(resRedefined)

  })

  raw_dat_res <- reactiveValues()
  
  observe({
    
    if(length(query()) != 0) {
    
    keyphrase <- as.character(query()[[1]])
    
    loadCaches(keyphrase, assignToVariable = "cipher", loadEnvir = globalenv(), cacheDir = "cache")
    
    cipher <- get("cipher", envir = globalenv())
    
    # keyphrase
    key <- hash(charToRaw(keyphrase))
    
    dat <- data_decrypt(cipher, key)
    
    raw_dat_res$raw_dat <- unserialize(dat)
    
    # disable all buttons in header when query is good
    shinyjs::disable("run")
    shinyjs::disable("userset")
    shinyjs::disable("universe")
    shinyjs::disable("loladb")
    shinyjs::disable("button_userset_example")
    shinyjs::disable("refgenome_core")
    shinyjs::disable("refgenome_ext")
    shinyjs::disable("defaultuniverse")
    shinyjs::disable("useruniverse")
    shinyjs::disable("checkbox")
    
    } else {
    
    raw_dat_res$raw_dat <- raw_dat_nocache()
    
    }

  })
  
  dat <- reactive({
    
      dat <- subset(raw_dat_res$raw_dat, maxRnk <= input$slider_rank_i)
      
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
    
    req(raw_dat_res$raw_dat)
    
    if(length(unique(raw_dat_res$raw_dat$userSet)) == 1) {
      
      unique(raw_dat_res$raw_dat$userSet)
      
    } else {
      
      c("All User Sets", unique(raw_dat_res$raw_dat$userSet))
      
    }
  }
  
  output$select_collection <- renderUI({
    
    req(raw_dat_res$raw_dat)
    
    list(
      HTML("<hr>"),
      selectInput("select_collection_i", 
                "Select Collection", 
                choices = c("All Collections", unique(raw_dat_res$raw_dat$collection)),
                selected = "All Collections"))
    
  })  
  
  output$slider_rank <- renderUI({
    
    req(raw_dat_res$raw_dat)
    
    sliderInput("slider_rank_i", 
                "Max Rank Cutoff", 
                min = min(raw_dat_res$raw_dat$maxRnk),
                max = max(raw_dat_res$raw_dat$maxRnk),
                value = quantile(raw_dat_res$raw_dat$maxRnk, probs = 20/nrow(raw_dat_res$raw_dat)))
    
  })  
  
  output$select_sort <- renderUI({
    
    req(raw_dat_res$raw_dat)
    
    selectInput("select_sort_i", 
                "Select Sort Column", 
                choices = names(raw_dat_res$raw_dat),
                selected = "meanRnk")
    
  })  
  
  output$select_userset <- renderUI({
    
    req(raw_dat_res$raw_dat)
    
    selectInput("select_userset_i", 
                "Select User Set", 
                choices = setchoices())
    
  })  
  
  # barplots
    
  # odds ratio
  
  output$slider_oddsratio <- renderUI({
    
    req(raw_dat_res$raw_dat)
    
    sliderInput("slider_oddsratio_i",
                "Odds Ratio Cutoff",
                min = round(min(raw_dat_res$raw_dat$oddsRatio), 3),
                max = round(max(raw_dat_res$raw_dat$oddsRatio), 3),
                value = round(min(raw_dat_res$raw_dat$oddsRatio), 3))
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
    
    req(raw_dat_res$raw_dat)
    
    sliderInput("slider_support_i",
                "Support Cutoff",
                min = round(min(raw_dat_res$raw_dat$support), 3),
                max = round(max(raw_dat_res$raw_dat$support), 3),
                value = round(min(raw_dat_res$raw_dat$support), 3))
    
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
    
    req(raw_dat_res$raw_dat)
    
    sliderInput("slider_pvalue_i", 
                "P Value Cutoff", 
                min = round(min(raw_dat_res$raw_dat$pValueLog), 3), 
                max = round(max(raw_dat_res$raw_dat$pValueLog), 3),
                value = round(min(raw_dat_res$raw_dat$pValueLog), 3))
    
    
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
