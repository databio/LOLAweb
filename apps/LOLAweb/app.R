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
library(shinyBS)

setCacheDir("cache")

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    # javascript for redirect to results view
    tags$script("Shiny.addCustomMessageHandler('redirect', 
                function(result_url) {window.location = result_url;});")   
  ),
  
  titlePanel(title = HTML("<img src='LOLA-logo.png' alt='LOLA logo' width='200'>"),
             windowTitle = "LOLA"),
  
      fluidRow(
        column(4,
               # need the bootstrap button invocation for popovers to work
               # but we don't actually want to show it
               shinyjs::hidden(bsButton("renderButton", "Render")),
               tags$div(
                 h3("#1 Select User Set",
                    # tags$a(href = "http://code.databio.org/LOLA/articles/gettingStarted.html", 
                    #        icon("question-circle-o"), 
                    #        target = "blank"))
                    actionLink("infouserset", "", icon = icon("question-circle-o")))
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
               ),
               tags$div(
               selectInput("refgenome", label = "Reference Genome", choices = c("hg19", "hg38", "mm10")),
               style = "margin-top:30px;")
        ),
        column(4,
               tags$div(
                 h3("#2 Select Universe",
                    # tags$a(href = "http://code.databio.org/LOLA/articles/choosingUniverse.html", 
                    #        icon("question-circle-o"), 
                    #        target = "blank"))
                 actionLink("infouniverse", "", icon = icon("question-circle-o")))
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
                    # tags$a(href = "http://databio.org/regiondb", 
                    #        icon("question-circle-o"), 
                    #        target = "blank"))
                 actionLink("infodb", "", icon = icon("question-circle-o")))
               ),
               # HTML(disabledbutton),
               uiOutput("loladbs")
               # conditionalPanel(condition = "input.loladb == 'Core'",
               #                  selectInput("refgenome_core", 
               #                              "Reference Genome", 
               #                              choices = c("hg19","hg38", "mm10"))),
               # conditionalPanel(condition = "input.loladb == 'Extended'",
               #                  selectInput("refgenome_ext", 
               #                              "Reference Genome", 
               #                              choices = c("hg19","hg38")))
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
      ),
  # footer text with SOMRC link
  tags$footer(HTML("<a href = 'https://somrc.virginia.edu'><i class='fa fa-bolt'></i> Powered By SOMRC</a>"), align = "right", style = " bottom:0; width:100%; height:10px; padding: 10px; z-index: 1000;"
  )
)

server <- function(input, output, session) {
  
  # popover text
  # user sets
  addPopover(session=session, 
             id="infouserset", 
             title="User Sets", 
             content="<p>The User Set is your set of genomic regions that you want to test for overlap with the database. Upload a file in <a href = 'https://genome.ucsc.edu/FAQ/FAQformat.html' target 'blank'>BED format</a> (really, it just needs the first 3 columns: chr, start, and end). You can also drag and drop multiple files and they will be analyzed simultaneously!<button type='button' id='close' class='close' onclick='$(&quot;#infouserset&quot;).popover(&quot;hide&quot;);'>&times;</button></p>", 
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  # universes
  addPopover(session=session, 
             id="infouniverse", 
             title="Universes", 
             content="<p>The universe is your background set of regions. You should think of the universe as the set of regions you tested for possible inclusion in your user sets; or, in other words, it is the restricted background set of regions that were tested, including everything in your regions of interest as well as those that did not get included. We recommend you upload a universe that is specific to the query set, but we also provide a few basic region sets (like tiling regions, or the set of all active DNaseI hypersensitive elements from the ENCODE project). The choice of universe can have a drastic affect on the results of the analysis, so it may also be worth running LOLA few times with different universe sets. For further information, there are more details in the <a href = 'http://code.databio.org/LOLA/articles/choosingUniverse.html' target='blank'>LOLA documentation</a>.<button type='button' id='close' class='close' onclick='$(&quot;#infouniverse&quot;).popover(&quot;hide&quot;);'>&times;</button></p>", 
             placement = "bottom",
             trigger = "click", 
             options = NULL)

  # region dbs
  addPopover(session=session, 
             id="infodb", 
             title="Region Databases", 
             # html for content with JS at the bottom to close popup
             content="<p>We have provided a few different general-purpose databases. We recommend starting with the Core database, but there are also a few other more specific options if you want to extend your analysis. Further details about what is contained in each database can be found in the <a href = 'http://databio.org/regiondb' target = 'blank'>documentation on LOLA region databases</a>. <button type='button' id='close' class='close' onclick='$(&quot;#infodb&quot;).popover(&quot;hide&quot;);'>&times;</button></p>", 
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  # addPopover(session=session, 
  #            id="infocollection", 
  #            title="collection", 
  #            # html for content with JS at the bottom to close popup
  #            content="<p>LOLA databases are made up of one or more sub-collections of region set. Using this drop-down, you can filter your plots and tables to show only the results from one of these collections at a time. <button type='button' id='close' class='close' onclick='$(&quot;#infocollection&quot;).popover(&quot;hide&quot;);'>&times;</button></p>", 
  #            placement = "bottom",
  #            trigger = "click", 
  #            options = NULL)
  
  
  # reactive values
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
      rawdat_nocache()
    },
    message = function(m) {
      shinyjs::html(id = "messages", html = m$message, add = FALSE)
        
    })
    
    shinyjs::hide("messages")
    shinyjs::hide("gear")
    
    showModal(modalDialog(
      title = "Results",
      HTML(paste0("You are about to be re-directed to your results:<br>",
      result_url()$result_url)
      )
    ))
    
    Sys.sleep(3)
    
    # initiate redirect
    session$sendCustomMessage(type = "redirect", result_url()$link)
    
  })
  
  result_url <- reactive({
    
    baseurl <- session$clientData$url_hostname
    port <- session$clientData$url_port
    pathname <- session$clientData$url_pathname
    
    # logic to remove phantom : when running at port 80
    if (port == 80) {
      
      link <- paste0("http://", baseurl, port, pathname, "?key=", keyphrase())
      
    } else {
      
      link <- paste0("http://", baseurl, ":", port, pathname, "?key=", keyphrase())
      
    }
  
    result_url <- paste0("<a href = '", link, "' target = 'blank'>", link, "</a>")
    
    # paste0("http://", baseurl, ":", port, pathname, "?key=", keyphrase())
    list(link = link,
         result_url = result_url)
    
  })
  
  output$link <- renderText({
    
    result_url()$result_url
    
  })
  
  
  output$loladbs <- renderUI({
    
    if(input$refgenome == "mm10") {
      
      selectInput("loladb", label = "", choices = c("Core"))
      
    } else {
      
      selectInput("loladb", label = "", choices = c("Core", "LOLAJaspar", "LOLARoadmap"))
      
    }
    
  })
  
  output$universe <- renderUI({

    if(input$checkbox) {
        
      fileInput("useruniverse", "Upload Universe",
                accept = c(".bed"))
      
      } else {
        
        selectInput("defaultuniverse", 
                    label = "Select Pre-Loaded Universe",
                    choices = list.files(paste0("universes/", input$refgenome)))
                    # choices = list.files("universes"))
  
        
      }
      
    })
    
  keyphrase <- eventReactive(input$run, {
      
    paste0(sample(c(LETTERS,1:9), 15), collapse = "")
      
    })
    
    rawdat_nocache <- eventReactive(input$run, {
    
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

        datapath <- paste0("universes/", input$refgenome, "/", input$defaultuniverse)

        userUniverse = read.table(file = datapath, header = F)

      }
      colnames(userUniverse) <- c('chr','start','end','id','score','strand')
      userUniverse <- with(userUniverse, GRanges(chr, IRanges(start+1, end), strand, score, id=id))

      userSetsRedefined =	redefineUserSets(userSets, userUniverse)

      # load region data for each reference genome
      dbPath = paste("reference",
                      input$loladb,
                      input$refgenome,
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

      simpleCache(cacheName = keyphrase(), 
                  instruction = { cipher },
                  noload = TRUE)

      return(resRedefined)

  })

  rawdat_res <- reactiveValues()
  
  observe({
    
    cachenames <- tools::file_path_sans_ext(list.files("cache"))
    
    if(length(query()) != 0 && query()[[1]] %in% cachenames) {
    
    keyphrase <- as.character(query()[[1]])
    
    env <- new.env()
    
    loadCaches(keyphrase, assignToVariable = "cipher", loadEnvir = env, cacheDir = "cache")
    
    cipher <- get("cipher", envir = env)
    
    # keyphrase
    key <- hash(charToRaw(keyphrase))
    
    dat <- data_decrypt(cipher, key)
    
    rawdat_res$rawdat <- unserialize(dat)
    
    # disable all buttons in header when query is good
    shinyjs::disable("run")
    shinyjs::disable("userset")
    shinyjs::disable("universe")
    shinyjs::disable("loladb")
    shinyjs::disable("button_userset_example")
    shinyjs::disable("refgenome")
    shinyjs::disable("defaultuniverse")
    shinyjs::disable("useruniverse")
    shinyjs::disable("checkbox")
    
    } else {
    
    rawdat_res$rawdat <- rawdat_nocache()
    
    }

  })
  
  dat <- reactive({
    
      dat <- subset(rawdat_res$rawdat, maxRnk <= input$slider_rank_i)
      
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
    
    req(rawdat_res$rawdat)
    
    if(length(unique(rawdat_res$rawdat$userSet)) == 1) {
      
      unique(rawdat_res$rawdat$userSet)
      
    } else {
      
      c("All User Sets", unique(rawdat_res$rawdat$userSet))
      
    }
  }
  
  output$select_collection <- renderUI({
    
    req(rawdat_res$rawdat)
    
    list(
      HTML("<hr>"),
      # tags$div(
        selectInput("select_collection_i", 
                    "Select Collection", 
                    choices = c("All Collections", unique(rawdat_res$rawdat$collection)),
                    selected = "All Collections")
           # actionLink("infocollection", "", icon = icon("question-circle-o")))
      # )
    )
  })  
  
  output$slider_rank <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_rank_i", 
                "Max Rank Cutoff", 
                min = min(rawdat_res$rawdat$maxRnk),
                max = max(rawdat_res$rawdat$maxRnk),
                value = quantile(rawdat_res$rawdat$maxRnk, probs = 20/nrow(rawdat_res$rawdat)))
    
  })  
  
  output$select_sort <- renderUI({
    
    req(rawdat_res$rawdat)
    
    selectInput("select_sort_i", 
                "Select Sort Column", 
                choices = names(rawdat_res$rawdat),
                selected = "meanRnk")
    
  })  
  
  output$select_userset <- renderUI({
    
    req(rawdat_res$rawdat)
    
    selectInput("select_userset_i", 
                "Select User Set", 
                choices = setchoices())
    
  })  
  
  # barplots
    
  # odds ratio
  
  output$slider_oddsratio <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_oddsratio_i",
                "Odds Ratio Cutoff",
                min = round(min(rawdat_res$rawdat$oddsRatio), 3),
                max = round(max(rawdat_res$rawdat$oddsRatio), 3),
                value = round(min(rawdat_res$rawdat$oddsRatio), 3))
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
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_support_i",
                "Support Cutoff",
                min = round(min(rawdat_res$rawdat$support), 3),
                max = round(max(rawdat_res$rawdat$support), 3),
                value = round(min(rawdat_res$rawdat$support), 3))
    
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
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_pvalue_i", 
                "P Value Cutoff", 
                min = round(min(rawdat_res$rawdat$pValueLog), 3), 
                max = round(max(rawdat_res$rawdat$pValueLog), 3),
                value = round(min(rawdat_res$rawdat$pValueLog), 3))
    
    
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
