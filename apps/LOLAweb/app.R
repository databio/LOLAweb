options(shiny.maxRequestSize=100*1024^2)

source("themes.R")
source("misc.R")

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)
library(DT)
library(simpleCache)
library(sodium)
library(shinyBS)
library(GenomicDistributions)

setCacheDir("cache")


ui <- list( 
  
  # # need empty fluid page at top with height 0 to preserve window title
  div(
    fluidPage(
      list(tags$head(HTML('<link rel="icon", href="LOLAweb-logo.png",
                          type="image/png" />'))),
      div(
        style="height:0px; padding:0px; width: '100%'",
        titlePanel(title="", windowTitle="LOLA")
        )
      ),
    style = "height:0px"),
  
  navbarPage(title = div(a(img(src = "LOLAweb-logo-cropped.png", style = "width:150px"), href = "/"), ""),
  
  tabPanel("Run",
  shinyjs::useShinyjs(),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    # javascript for redirect to results view
    tags$script("Shiny.addCustomMessageHandler('redirect', 
                function(result_url) {window.location = result_url;});"),
    tags$link(rel="shortcut icon", href="favicon.ico")

  ),
      fluidRow(
        shinyjs::hidden(div(HTML("<div class='alert alert-warning'>
          All <strong>Run</strong> inputs are disabled while LOLAweb has results loaded.<br>Navigate to <strong>Results</strong> to view current result output.<br>To execute a new run, <a href= '/'>reload</a> LOLAweb.</div>"), id = "noinputmsg"))
          
      ),
      fluidRow(
        column(4,
               # need the bootstrap button invocation for popovers to work
               # but we don't actually want to show it
               shinyjs::hidden(bsButton("renderButton", "Render")),
               tags$div(
                h3("#1 Upload regions of interest",
                    actionLink("infouserset", "", icon = icon("question-circle-o")))
               ),
               selectInput("refgenome", label = "Reference genome", choices = c("hg19", "hg38", "mm10")),
               shinyjs::hidden(selectInput("defaultuserset", 
                                           label = "Select pre-loaded region set",
                                           choices = list.files("userSets"))),
               fileInput("userset", "Choose BED file(s)",
                         multiple = TRUE,
                         accept = c(".bed")),
                tags$div(actionButton("button_userset_example",
                            "Load example data"), style="margin-top:-15px"),
               shinyjs::hidden(
                 actionButton("button_userset_upload",
                            "Browse for data to upload")
               )
        ),
        column(4,
               tags$div(
                 h3("#2 Select background universe",
                 actionLink("infouniverse", "", icon = icon("question-circle-o")))
               ),
               uiOutput("universe"),
               checkboxInput("checkbox", 
                             label = "Check here to upload your own universe",
                             value = FALSE)
               ),
        column(4, 
               tags$div(
                 h3("#3 Select region database",
                 actionLink("infodb", "", icon = icon("question-circle-o")))
               ),
               uiOutput("loladbs"),
               actionButton("run",
                            "RUN LOLA", 
                            class = "runLOLA")
        ),
      class = "headerBox"),
  fluidRow(
    column(2,
           htmlOutput("gear")),
    column(10,
           htmlOutput("messages"))
  )),
  tabPanel("Results",
        fluidRow(div(HTML("<div class='alert alert-warning'><strong>Results</strong> is currently empty.<br>To generate result output, visit <strong>Run</strong> or view <a href = '?key=CSD6HEAN2RJUOVK'>sample results</a>.</div>"), id = "noresultsmsg")
             ),
        fluidRow(
        column(2,
               shinyjs::hidden(
                 tags$div(
                   h4("Display Options",
                      actionLink("infodisplay", "", icon = icon("question-circle-o"))),
                   id = "infodisplay_div"
                 )
               ),
               shinyjs::hidden(htmlOutput("gear2")),
               uiOutput("slider_rank"),
               uiOutput("slider_oddsratio"),
               uiOutput("slider_support"),
               uiOutput("slider_pvalue"),
               uiOutput("select_collection"),
               uiOutput("select_sort"),
               uiOutput("select_userset")
          ),
        column(10,
               tags$h4(htmlOutput("link")),
               shinyjs::hidden(
                 tags$div(
                   h2("LOLA Results",
                      actionLink("infoplot", "", icon = icon("question-circle-o"))),
                   id = "infoplot_div")),
               shinyjs::hidden(htmlOutput("messages2"))
               ),
        column(5,
               conditionalPanel(condition = "output.res",
                                h4("Odds Ratio"),
                                downloadButton("oddsratio_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("oddsratio_plot"),
               conditionalPanel(condition = "output.res",
                                h4("Support"),
                                downloadButton("support_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("support_plot")
        ),
        column(5,
               conditionalPanel(condition = "output.res",
                                h4("P Value"),
                                downloadButton("pvalue_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("pvalue_plot"),
               conditionalPanel(condition = "output.res",
                                h4("Distribution over genome"),
                                downloadButton("distrib_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("distrib_plot"),
               conditionalPanel(condition = "output.res",
                                h4("Distance to TSS"),
                                downloadButton("dist_plot_dl",
                                               label = "Download Plot",
                                               class = "dt-button")),
               plotOutput("dist_plot")
        )
        ),
      fluidRow(
        column(DT::dataTableOutput("res"), width = 12) 
      )
  ),
  tabPanel("About",
    includeHTML("about.html")
  ),
  footer = div(HTML("<div>Powered by <a href = 'https://somrc.virginia.edu' target ='blank'>SOMRC</a><br>A Project of the <a href ='http://databio.org/' target = 'blank'>Sheffield Lab</a><br>Source code on <a href ='https://github.com/databio/LOLAweb' target = 'blank'>Github</a><br>Try it yourself with <a href='https://github.com/databio/LOLAweb/blob/master/docker/README.md' target = 'blank'>Docker</a></div>"), align = "right", style = " bottom:0; width:100%; height:10px; padding: 10px; padding-bottom:20px; z-index: 1000;"),
  id = "mainmenu"
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
  
  # collection
  addPopover(session=session,
             id="infocollection",
             title="Collections",
             # html for content with JS at the bottom to close popup
             content="<p>LOLA databases are made up of one or more sub-collections of region set. Using this drop-down, you can filter your plots and tables to show only the results from one of these collections at a time. <button type='button' id='close' class='close' onclick='$(&quot;#infocollection&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click",
             options = NULL)
  
  # LOLA results
  addPopover(session=session,
             id="infoplot",
             title="LOLA Results",
             # html for content with JS at the bottom to close popup
             content="<p>These barplots show the highest-ranking region sets from the database. The higher scores indicate more overlap with your query set. The results are scored using 3 statistics: Support is the raw number of regions that overlapped between your query set and the database set. LogPVal and LogOdds are the results of a Fisher's exact test scoring the significance of that overlap.</p><p>We rank each region set from the database for each of these 3 scores, and you can see the raw scores and the ranks in the table below. You can also see the maximum and mean ranks across all 3 categories. <button type='button' id='close' class='close' onclick='$(&quot;#infoplot&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click",
             options = NULL)
  
  # display slider options
  addPopover(session=session,
             id="infodisplay",
             title="Display Options",
             # html for content with JS at the bottom to close popup
             content="<p>These sliders can be used to adjust the number of results displayed on the plots. You may adjust cutoffs for any of the 3 statistics, as well as for the combined maximum rank, which prioritizes region sets that score well in all 3 categories.</p><p>LOLA databases are made up of one or more sub-collections of region set. Using this drop-down, you can filter your plots and tables to show only the results from one of these collections at a time. <button type='button' id='close' class='close' onclick='$(&quot;#infodisplay&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click",
             options = NULL)
  
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
        
      fileInput("useruniverse", "Choose universe file",
                accept = c(".bed"))
      
      } else {
        
        selectInput("defaultuniverse", 
                    label = "Select pre-loaded universe",
                    choices = list.files(paste0("universes/", input$refgenome)))
      
      }
      
    })
    
  keyphrase <- eventReactive(input$run, {
      
    paste0(sample(c(LETTERS,1:9), 15), collapse = "")
      
    })
    
    rawdat_nocache <- eventReactive(input$run, {
    
      message("Calculating region set enrichments ...")
      userSets <- list()

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
      
      # garbage collection
      gc()

      resRedefined = runLOLA(userSetsRedefined,
                             userUniverse,
                             regionDB,
                             cores=4)

      # need to make sure user set is discrete even if coded as number
      resRedefined$userSet = as.character(resRedefined$userSet)
      
      # calculate distribution over chromosomes for plotting
      genDist = aggregateOverGenomeBins(userSets, input$refgenome)

      # calculate distances to TSSs
      TSSDist = TSSDistance(userSets, input$refgenome)

      # create named list of multiple objects for plotting
      res = list(resRedefined = resRedefined,
                 genDist = genDist,
                 TSSDist = TSSDist)
      
      # caching
      # need to call keyphrase from reactive above because it is used to construct link
      key <- hash(charToRaw(keyphrase()))
      msg <- serialize(res,
                       connection = NULL)

      cipher <- data_encrypt(msg, key)

      simpleCache(cacheName = keyphrase(), 
                  instruction = { cipher },
                  noload = TRUE)

      return(list(resRedefined = resRedefined,
                  genDist = genDist,
                  TSSDist = TSSDist))

  })

  rawdat_res <- reactiveValues()
  
  observe({
    
    cachenames <- tools::file_path_sans_ext(list.files("cache"))
    
    if(length(query()) != 0 && !query()[[1]] %in% cachenames) {
      
      showModal(modalDialog(
        title = "Bad Request",
        HTML(paste0("The cache '",
                     query()[[1]], 
                    "' does not exist."))
        )
      )
    }
    
  })
  
  plot_render <- reactiveValues(state = FALSE)
  
  observe({
    
    cachenames <- tools::file_path_sans_ext(list.files("cache"))
    
    if(length(query()) != 0 && query()[[1]] %in% cachenames) {
    
    keyphrase <- as.character(query()[[1]])
    
    env <- new.env()
    
    simpleCache(keyphrase, assignToVariable = "cipher", loadEnvir = environment(), cacheDir = "cache")
    
    cipher <- get("cipher", envir = env)
    
    # keyphrase
    key <- hash(charToRaw(keyphrase))
    
    dat <- data_decrypt(cipher, key)
    
    res <- unserialize(dat)
    
    rawdat_res$rawdat <- res$resRedefined
    
    rawdat_res$genDist <- res$genDist

    rawdat_res$TSSDist <- res$TSSDist
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
    
    # show results message on run tab
    shinyjs::show("noinputmsg")
    
    # hide no results message
    shinyjs::hide("noresultsmsg")
    
    updateNavbarPage(session, "mainmenu",
                      selected = "Results")
    
    shinyjs::show("gear2")
    shinyjs::show("messages2")
    
    # show help text for results sliders and plots
    shinyjs::show("infoplot_div")
    shinyjs::show("infodisplay_div")
    
    } else {
    
    rawdat_res$rawdat <- rawdat_nocache()$resRedefined
    rawdat_res$genDist <- rawdat_nocache()$genDist
    
    }

  })
  
  observe({
    
    if(!plot_render$state) {
      
      shinyjs::html(id = "gear2", html = "<i class='fa fa-4x fa-spin fa-cog'></i>", add = FALSE)
      shinyjs::html(id ="messages2", html = "Rendering plots ...", add = FALSE)
      
    } else {
      
      shinyjs::hide("gear2")
      shinyjs::hide("messages2")
      
    }
    
  })

  dat <- reactive({
    
      dat <- subset(rawdat_res$rawdat, maxRnk <= input$slider_rank_i)
      
      dat <- subset(dat,
                    oddsRatio >= input$slider_oddsratio_i &
                      support >= input$slider_support_i &
                      pValueLog >= input$slider_pvalue_i)
      
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
    
        selectInput("select_collection_i", 
                    "Select Collection", 
                    choices = c("All Collections", unique(rawdat_res$rawdat$collection)),
                    selected = "All Collections")
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
    
    sortcols <- c("pValueLog", "oddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk")
    
    selectInput("select_sort_i", 
                "Select Sort Column", 
                choices = sortcols,
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
  
  # call plot
  output$oddsratio_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    plot_render$state <- TRUE
    
    plot_input(dat(), "oddsRatio", "Odds Ratio", input$select_sort_i)

  })
  
  # download handler
  output$oddsratio_plot_dl <- downloadHandler(
    filename = function() { paste("oddsratio", 
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = plot_input(dat(), "oddsRatio", "Odds Ratio", input$select_sort_i)
, device = "pdf")
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
  
  output$support_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    plot_input(dat(), "support", "Support", input$select_sort_i)
    
  })
  
  # download handler
  output$support_plot_dl <- downloadHandler(
    filename = function() { paste("support", 
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = plot_input(dat(), "support", "Support", input$select_sort_i)
, device = "pdf")
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
  
  output$pvalue_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    plot_input(dat(), "pValueLog", "P Value (log scale)", input$select_sort_i)
    
    
  })
  
  output$pvalue_plot_dl <- downloadHandler(
    filename = function() { paste("pvalue",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = plot_input(dat(), "pValueLog", "P Value (log scale)", input$select_sort_i)
, device = "pdf")
    }
  )
  
  # genomic distribution plot
  # set up function
  distrib_plot_input <- function() {
    
    genDist <- rawdat_res$genDist
    
    plotGenomeAggregate(genomeAggregate = genDist)
    
  }
  
  output$distrib_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    distrib_plot_input()
    
  })
  
  # feature distance plot
  dist_plot_input <- function() {
    
    TSSDist <- rawdat_res$TSSDist
    
    plotFeatureDist(TSSDist, featureName="TSS")
    
  }
  
  output$dist_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    dist_plot_input()
    
  })

  # download handler
  output$distrib_plot_dl <- downloadHandler(
    filename = function() { paste("gendist",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = distrib_plot_input(), device = "pdf")
    }
  )
  
  output$dist_plot_dl <- downloadHandler(
    filename = function() { paste("tssdist",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = dist_plot_input(), device = "pdf")
    }
  )
  


  # data table
  output$res <- DT::renderDataTable({
    
    req(input$select_sort_i)
    
    if(grepl("rnk", input$select_sort_i, ignore.case = TRUE)) {
    
      dat <- dat()[order(eval(parse(text = input$select_sort_i)), decreasing = FALSE)]
      
    } else {
      
      dat <- dat()[order(eval(parse(text = input$select_sort_i)), decreasing = TRUE)]
      
    }

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