options(shiny.maxRequestSize=100*1024^2)

source("themes.R")
source("misc.R")
source("disabler.R")

db_url <- "mongodb://lw-db:27017"
db_name <- "lolaweb"

# where do the caches live?
# cache_dir <- "cache"

# establish connection to jobdb
con <- shinyqueue::connect(db_url = db_url, db_name = db_name)

library(shiny)
library(LOLA)
library(ggplot2)
library(GenomicRanges)
library(DT)
library(simpleCache)
library(sodium)
library(shinyBS)
library(GenomicDistributions)
library(plotly)

setCacheDir(cacheDir)

# cache_dir <- cacheDir
  
# get lolaweb version
lw_version <- system(command = "git rev-parse HEAD | cut -c1-9", intern = TRUE)

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
    # tags$script('Shiny.addCustomMessageHandler("jsCode",
    #                               function(message) {
    #                                 console.log(message)
    #                                 eval(message.code);
    #                               }
    # );
    # '),
    tags$link(rel="shortcut icon", href="favicon.ico")

  ),
      fluidRow(
        shinyjs::hidden(div(HTML("<div class='alert alert-warning'>
          All <strong>Run</strong> inputs are disabled while LOLAweb has results
          loaded.<br>Navigate to <strong>Results</strong> to view current result
          output.<br>To execute a new run, <a href= '/'>reload</a>
          LOLAweb.</div>"), id = "noinputmsg"))
          
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
               uiOutput("refgenome"),
               shinyjs::hidden(
                 uiOutput("defaultuserset")),
               fileInput("userset", "Choose BED file(s)",
                         multiple = TRUE,
                         accept = c(".bed")),
                tags$div(actionButton("button_userset_example",
                            "Load example data"), style="margin-top:-15px"),
               shinyjs::hidden(
                 actionButton("button_userset_upload",
                            "Browse for data to upload")
               ),
               textOutput("testn")
        ),
        column(4,
               tags$div(
                 h3("#2 Select background universe",
                 actionLink("infouniverse", "", icon = icon("question-circle-o")))
               ),
               uiOutput("universe_opts"),
               uiOutput("universe")
               ),
        column(4, 
               tags$div(
                 h3("#3 Select region database",
                 actionLink("infodb", "", icon = icon("question-circle-o")))
               ),
               uiOutput("loladbs"),
               actionButton("run",
                            "RUN LOLA", 
                            class = "runLOLA"),
               HTML("<div id='samplereslink' style='padding-top:15px; padding-left:5px;'><a href = '?key=M2LZOJXEQSKG98A'>Sample Results</a></div>")
        ),
  id = "runInputs"),
  
  fluidRow(
    column(2,
           htmlOutput("gear")),
    column(10,
           htmlOutput("messages"))
  )),
  tabPanel("Results",
           fluidRow(div(HTML("<div class='alert alert-warning'><strong>Results</strong> is currently empty.<br>To generate result output, visit <strong>Run</strong> or view <a href = '?key=M2LZOJXEQSKG98AE'>sample results</a>.</div>"), id = "noresultsmsg")
           ),
           fluidRow(
             column(10,
                    tags$h4(htmlOutput("link")),
                    shinyjs::hidden(htmlOutput("gear2"))
             )
           ),
           fluidRow(column(2,
                    shinyjs::hidden(
                      tags$div(
                        h4("Display Options",
                           actionLink("infodisplay", "", icon = icon("question-circle-o"))),
                        id = "infodisplay_div"
                             )),
                    uiOutput("slider_rank"),
                    uiOutput("slider_oddsratio"),
                    uiOutput("slider_support"),
                    uiOutput("slider_pvalue"),
                    uiOutput("select_collection"),
                    uiOutput("select_sort"),
                    uiOutput("select_userset"),
                    shinyjs::hidden(
                      downloadButton("all_plots_dl",
                                     label = "Download all plots",
                                     class = "dt-button"))),
                    column(10,
                           shinyjs::hidden(
                             div(
                             tabsetPanel(type = "tabs",
                                       tabPanel("Scatterplot",
                                                shinyjs::hidden(
                                                  div(
                                                    h4("Scatterplot", class ="plot_header"),
                                                    downloadButton("scatterplot_dl",
                                                                   label = "PDF",
                                                                   class = "dt-button"),
                                                    id = "scatterhead"
                                                )),
                                                plotlyOutput("scatter")),
                                       tabPanel("Barplots",
        column(5,
                                h4("Odds ratio", class = "plot_header"),
                                downloadButton("oddsratio_plot_dl",
                                               label = "PDF",
                                               class = "dt-button"),
               plotOutput("oddsratio_plot"),
                                h4("Support (overlap count)", class = "plot_header"),
                                downloadButton("support_plot_dl",
                                               label = "PDF",
                                               class = "dt-button"),
               plotOutput("support_plot")
        
        ),
        column(5,
                                h4("P Value", class = "plot_header"),
                                downloadButton("pvalue_plot_dl",
                                               label = "PDF",
                                               class = "dt-button"),
               plotOutput("pvalue_plot"), HTML("<i>These barplots depict the
               top-ranking database region sets. The sort order for all 3 plots
               is determined by the sort column selected in the Display Options
               panel. At most, 50 region sets will be displayed. Entries with
               multiple small bars represent multiple database replicates with
               the same descriptive label. Further details on each comparison
               can be found in the Table tab. </i>")
        )
           ),
         
                                       
      tabPanel("Genomic distribution",
          HTML("<i>These plots come from the 
            <a href='http://code.databio.org/GenomicDistributions/'>GenomicDistributions</a> 
            package.</i>"),
               fluidRow(
                 column(10,
                        h4("Distribution over genome", class = "plot_header"),
                        downloadButton("distrib_plot_dl",
                                       label = "PDF",
                                       class = "dt-button"),
                        plotOutput("distrib_plot"))
               ),
               fluidRow(
                 column(5,
                        h4("Distance to TSS", class = "plot_header"),
                        downloadButton("dist_plot_dl",
                                       label ="PDF",
                                       class = "dt-button"),
                        plotOutput("dist_plot")),
                 column(5,
                        h4("Genomic partitions", class = "plot_header"),
                        downloadButton("part_plot_dl",
                                       label = "PDF",
                                       class = "dt-button"),
                        plotOutput("part_plot"))
               )
               ),
      tabPanel("Table",
               column(DT::dataTableOutput("res"), width = 12)
               ),
      tabPanel("Run summary",
                                h4("Run summary"),
                                tableOutput("run_sum"),
               style = "font-size:18px;")
           ),
      id = "result-tabs")))
  )
  ),
  tabPanel("About",
    includeMarkdown("about.md")
  ),
  footer = div(HTML(
    paste0("<div>
      Built by the <a href ='http://databio.org/'
    target = 'blank'>Sheffield Computational Biology Lab</a> and <a href = 'https://somrc.virginia.edu'
    target='blank'>SOMRC</a> at UVA. <br>View <a href
    ='https://github.com/databio/LOLAweb' target = 'blank'>source code on GitHub</a> or run it locally with <a
    href='https://github.com/databio/LOLAweb/blob/master/docker/README.md'
    target = 'blank'>our docker image</a>", "<br>LOLAweb version: <a href
    ='https://github.com/databio/lolaweb/commit/", lw_version, "'>", lw_version,
    "</a></div>")
    ), 
    align = "center", style = " bottom:0; width:100%; height:10px; padding: 10px; padding-bottom:20px; z-index: 1000;"),
  id = "mainmenu"
)
)

server <- function(input, output, session) {
  
  # popover text
  # user sets
  addPopover(session=session, 
             id="infouserset", 
             title="User query sets", 
             content="<p>The regions of interest, also called the <i>user query
             set</i>, is your set of genomic regions that you want to test for
             overlap with the database. Upload a file in <a href =
             'https://genome.ucsc.edu/FAQ/FAQformat.html' target 'blank'>BED
             format</a> (really, it just needs the first 3 columns: chr, start,
             and end). You can also drag and drop multiple files and they will
             be analyzed simultaneously!</p>
             
             <p>You can also load example data for each reference genome.
             Details about what regions are in these examples are on the About
             tab.

             <button type='button' id='close' class='close' onclick='$(&quot;#infouserset&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  # universes
  addPopover(session=session, 
             id="infouniverse", 
             title="Background universe", 
             content="<p>The universe set is your background set of regions. You
             should think of the universe as the set of regions you tested for
             possible inclusion in your user sets. You have 3 options: 1) we
             provide a few pre-loaded region sets (like tiling regions, or the
             set of all active DNaseI hypersensitive elements from the ENCODE
             project) that can be useful for many analyses. 2), you can upload
             your own universe, which lets you tailor the analysis to the query
             set; 3) if you upload more than 1 query set, we can build a
             universe by merging all of the uploaded query sets. The choice of
             universe can have a drastic affect on the results of the analysis,
             it may be worth running LOLA few times with different universe
             sets. For further information, there are more details on the
             <b>About tab</b>. <button type='button' id='close' class='close'
             onclick='$(&quot;#infouniverse&quot;).popover(&quot;h
             ide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)

  # region dbs
  addPopover(session=session, 
             id="infodb", 
             title="Region databases", 
             # html for content with JS at the bottom to close popup
             content="<p>We have provided a few different general-purpose
             databases. We recommend starting with the Core database, which
             contains data from the ENCODE project, the cistrome project, and
             other public data sources. We also provide databases containing
             Roadmap Epigenomics data and motif searches using motifs from the
             JASPAR motif database. Further details about what is contained in
             each database can be found on the About tab. <button type='button'
             id='close' class='close' oncli
             ck='$(&quot;#infodb&quot;).popover(&quot;hide&quot;);'>&times;</but
             ton></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  # collection
  addPopover(session=session,
             id="infocollection",
             title="Collections",
             # html for content with JS at the bottom to close popup
             content="<p>LOLA databases are made up of one or more sub-
             collections of region set. Using this drop-down, you can filter
             your plots and tables to show only the results from one of these
             collections at a time. <button type='button' id='close' class='close' onclick='$(&quot;#infocollection&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click",
             options = NULL)
  
  # LOLA results
  addPopover(session=session,
             id="infoplot",
             title="LOLA Results",
             # html for content with JS at the bottom to close popup
             content="<p>These barplots show the highest-ranking region sets
             from the database. The higher scores indicate more overlap with
             your query set. The results are scored using 3 statistics: Support
             is the raw number of regions that overlapped between your query set
             and the database set. LogPVal and LogOdds are the results of a
             Fisher's exact test scoring the significance of that
             overlap.</p><p>We rank each region set from the database for each
             of these 3 scores, and you can see the raw scores and the ranks in
             the table below. You can also see the maximum and mean ranks across
             all 3 categories.
             <button type='button' id='close' class='close' onclick='$(&quot;#infoplot&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
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
  exampleuserset <- reactiveValues(toggle = FALSE)
  
  # get url parameter for retrieval query key
  
  # query <- reactive({
  #   
  #   parseQueryString(session$clientData$url_search)
  #   
  #   })
  
  query <- shinyqueue::parse_query()
  
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
    
    # disable runLOLA button while processing
    shinyjs::runjs("$('#runInputs').addClass('disabledinputs');")
    
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
  
  output$refgenome <- renderUI({
      
      selectInput("refgenome", label = "Reference genome", choices = c("hg19", "hg38", "mm10", "mm9"))
      
  })
  
  output$loladbs <- renderUI({

    fl = Sys.glob(paste0(dbDir, "/*/", input$refgenome))
    dbs = gsub(paste0(dbDir, "/(.*)/", input$refgenome), "\\1", fl)

    selectInput("loladb", label = "", choices = dbs)
    
  })
  
  output$defaultuserset <- renderUI({
    
    selectInput("defaultuserset",
                label = "Select pre-loaded region set",
                choices = list.files(paste0(exampleDir, "/", input$refgenome)))
    
  })

  output$universe_opts <- renderUI({
    
    if(nfiles$n > 1) {
      
      radioButtons("universe_opts", "Universe", 
                   choiceValues = c("default", "build", "user"), 
                   choiceNames = c("Use pre-loaded universe", "Build universe with user sets", "Upload universe"))

    } else {
      
      HTML(disabledbutton)
      
    }
    
  })
  
  output$universe <- renderUI({
    
    req(input$universe_opts)

    if(input$universe_opts == "user") {
        
      fileInput("useruniverse", "Choose universe file",
                accept = c(".bed"))
      
      } else if (input$universe_opts == "default") {
        
        selectInput("defaultuniverse", 
                    label = "Select pre-loaded universe",
                    choices = list.files(paste0(universeDir, input$refgenome)))
      
      }
      
    })
    
  keyphrase <- eventReactive(input$run, {
      
    # paste0(sample(c(LETTERS,1:9), 15), collapse = "")
    shinyqueue::job_hash()
      
    })
    
    rawdat_nocache <- eventReactive(input$run, {
      # 
      # message("Calculating region set enrichments ...")
      # 
      # start_time <- Sys.time()
      # 
      # run_time <- system.time({
      #     
      #   # user set(s)
      #   userSets <- list()
      # 
      #   if(!exampleuserset$toggle) {
      # 
      #     message("Loading uploaded data")
      # 
      #     for (i in 1:length(input$userset[,1])) {
      # 
      #       userSet = readBed(input$userset[[i, 'datapath']])
      #       userSets[[i]] <- userSet
      # 
      #     }
      # 
      #     userSets = GRangesList(userSets)
      # 
      #     names(userSets) = input$userset[,"name"]
      # 
      #   } else {
      # 
      #     message("Loading example data")
      # 
      #     datapath <- paste0(exampleDir, input$refgenome, "/", input$defaultuserset)
      # 
      #     userSet = readBed(datapath)
      #     userSets[[1]] <- userSet
      # 
      #     userSets = GRangesList(userSets)
      # 
      #     names(userSets) = input$defaultuserset
      # 
      #   }
      #   
      #   # universe
      #   
      #   if(input$universe_opts == "user") {
      #     
      #     userUniverse <- readBed(file = input$useruniverse$datapath)
      #     
      #     universename <- input$useruniverse$datapath
      #     
      #   } else if (input$universe_opts == "default") {
      #     
      #     datapath <- paste0(universeDir, input$refgenome, "/", input$defaultuniverse)
      #     
      #     userUniverse <- readBed(file = datapath)
      #   
      #     universename <- input$defaultuniverse
      #     
      #   } else if (input$universe_opts == "build") {
      #     
      #     userUniverse <- buildRestrictedUniverse(userSets)
      #     
      #     universename <- paste0(c(names(userSets), "universe"), collapse = "_")
      #     
      #   }
      #   
      #   userSetsRedefined =	redefineUserSets(userSets, userUniverse)
      # 
      #   # set up path to databases for given reference genome
      #   dbPath = paste(dbDir,
      #                   input$loladb,
      #                   input$refgenome,
      #                  sep = "/"
      #                   )
      #   
      #   # load region data for the given reference genome
      #   regionDB = loadRegionDB(dbLocation=dbPath)
      #   
      #   # garbage collection
      #   gc()
      # 
      #   resRedefined = runLOLA(userSetsRedefined,
      #                          userUniverse,
      #                          regionDB,
      #                          cores=4)
      #   
      #   # need to make sure user set is discrete even if coded as number
      #   resRedefined$userSet = as.character(resRedefined$userSet)
      #   
      #   gd_data <- utils::data(package="GenomicDistributions")$results[,"Item"]
      #   
      #   # calculate distribution over chromosomes for plotting
      #   if (paste0("chromSizes_", input$refgenome) %in% gd_data) {
      #     
      #     genDist = aggregateOverGenomeBins(userSets, input$refgenome)
      #     
      #   } else {
      #     
      #     genDist = NULL
      #     
      #   }
      # 
      #   # calculate distances to TSSs
      #   if (paste0("TSS_", input$refgenome) %in% gd_data) {
      #     
      #     TSSDist = TSSDistance(userSets, input$refgenome)
      #   
      #     
      #   } else {
      #     
      #     TSSDist = NULL
      #     
      #   }
      #   
      #   # distribution of overlaps for a query set to genomic partitions
      #   if (paste0("geneModels_", input$refgenome) %in% gd_data) {
      #  
      #     gp = genomicPartitions(userSets, input$refgenome)
      #   
      #   } else {
      #     
      #     gp = NULL
      #     
      #   }
      #     
      # })
      # 
      # run_sum <- 
      #   data.frame(
      #     start_time = as.character(start_time),
      #     end_time = as.character(Sys.time()),
      #     run_time = paste0(round(run_time[3])+1, " seconds"),
      #     cache_name = keyphrase(),
      #     query_set = paste(gsub(".bed","",unique(resRedefined$userSet)), collapse = "\n"),
      #     genome = input$refgenome,
      #     universe = gsub(".bed", "", universename),
      #     region_db = input$loladb,
      #     commit = lw_version
      #     )
      # 
      #   # create named list of multiple objects for plotting
      # res = list(resRedefined = resRedefined,
      #            genDist = genDist,
      #            TSSDist = TSSDist,
      #            run_sum = run_sum,
      #            gp = gp)
      # 
      # # caching
      # # need to call keyphrase from reactive above because it is used to construct link
      # key <- hash(charToRaw(keyphrase()))
      # msg <- serialize(res,connection = NULL)
      # 
      # cipher <- data_encrypt(msg, key)
      # 
      # simpleCache(cacheName = keyphrase(), 
      #             instruction = { cipher },
      #             noload = TRUE)
      # 
      # return(res)
      
      for (i in 1:length(input$userset[,1])) {
        
        # # get basename of input file
        # bn <- basename(input$userset[[i, 'datapath']])
        # 
        # # overwrite with tempdir
        # input$userset[[i, 'datapath']]  <- paste0("~/VP/tempdir/", bn)
        
        file.copy(input$userset[[i, 'datapath']],
                   file.path(tempDir, input$userset[[i, 'name']]) )
        
        }
      
      
      shinyqueue::submit(con, 
                         job_type = "lolaweb",
                         cache_dir = cacheDir,
                         input = input,
                         job_id = keyphrase(),
                         encrypt = TRUE)

  })

  rawdat_res <- reactiveValues()
  
  observe({
    
    cachenames <- tools::file_path_sans_ext(list.files(cacheDir))
    
    if(length(query() != 0)) {
      
      processed <- query()[[1]] %in% cachenames
      submitted <- query()[[1]] %in% con$find()$job_id
      bad_query <- length(query()) != 0 && !processed && !submitted
    
    if(bad_query) {
      # message("Cachenames: ", cachenames)
      # message("cacheDir: ", cacheDir)
      showModal(modalDialog(
        title = "Bad Request",
        HTML(paste0("The cache '",
                     query()[[1]], 
                    "' does not exist."))
        )
      )
    } else if (submitted && !processed) {
      
      showModal(modalDialog(
        title = "In progess",
        HTML(paste0("The cache '",
                    query()[[1]], 
                    "' does not exist yet.<br>The job is currently in the queue to be processed."))
      )
      )
    }
    }
    
  })
  
  # capture the number of user sets that are being uploaded
  # use this to toggle radio buttons for building a universe 
  nfiles <- reactiveValues(n = 0)
  
  observe({
    
    nfiles$n <- length(input$userset[,1])
    
  })
  
  plot_render <- reactiveValues(state = FALSE)
  
  observe({
    
    cachenames <- tools::file_path_sans_ext(list.files(cacheDir))
    
    if(length(query()) != 0 && query()[[1]] %in% cachenames) {
    
    keyphrase <- as.character(query()[[1]])
    
    env <- new.env()
    
    simpleCache(keyphrase, assignToVariable = "cipher", loadEnvir = environment(), cacheDir=cacheDir)
    
    cipher <- get("cipher", envir = env)
    
    # keyphrase
    key <- hash(charToRaw(keyphrase))
    
    dat <- data_decrypt(cipher, key)
    
    res <- unserialize(dat)
    
    rawdat_res$rawdat <- res$resRedefined
    
    rawdat_res$genDist <- res$genDist

    rawdat_res$TSSDist <- res$TSSDist
    
    rawdat_res$gp  <- res$gp
    
    rawdat_res$run_sum <- res$run_sum
    
    # disable all buttons in header when query is good
    shinyjs::runjs("$('#runInputs').addClass('disabledinputs');")
    
    # show results message on run tab
    shinyjs::show("noinputmsg")
    
    # hide no results message
    shinyjs::hide("noresultsmsg")
    
    updateNavbarPage(session, "mainmenu",
                      selected = "Results")
    
    shinyjs::show("gear2")
    shinyjs::show("result-tabs")
    
    # show help text for results sliders and plots
    shinyjs::show("infoplot_div")
    shinyjs::show("infodisplay_div")
    shinyjs::show("run_sum")
    
    
  
    } else {
    
    rawdat_res$rawdat <- rawdat_nocache()$resRedefined
    rawdat_res$genDist <- rawdat_nocache()$genDist
    rawdat_res$TSSDist <- rawdat_nocache()$TSSDist
    rawdat_res$gp <- rawdat_nocache()$gp
    rawdat_res$run_sum <- rawdat_nocache()$run_sum
    
    }

  })
  
  observe({
    
    if(!plot_render$state) {
      
      shinyjs::html(id = "gear2", html = "<i class='fa fa-4x fa-spin fa-cog'></i>", add = FALSE)

    } else {
      
      shinyjs::hide("gear2")
      shinyjs::show("scatterhead")
      shinyjs::show("all_plots_dl")
      
    
    }
    
  })

  dat <- reactive({
    
      dat <- subset(rawdat_res$rawdat, maxRnk <= input$slider_rank_i)
      
      dat <- subset(dat,
                    oddsRatio >= input$slider_oddsratio_i &
                      support >= input$slider_support_i &
                      pValueLog >= input$slider_pvalue_i)
      
      if (input$select_collection_i != "All collections") {
        
         dat <- subset(dat, collection == input$select_collection_i)
         
      }
      
      if (input$select_userset_i != "All user sets") {
  
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
      
      c("All user sets", unique(rawdat_res$rawdat$userSet))
      
    }
  }
  
  output$select_collection <- renderUI({
    
    req(rawdat_res$rawdat)
    
        selectInput("select_collection_i", 
                    "Filter by collection", 
                    choices = c("All collections", unique(rawdat_res$rawdat$collection)),
                    selected = "All collections")
  })  
  
  output$slider_rank <- renderUI({
    
    req(rawdat_res$rawdat)
    sliderInput("slider_rank_i", 
                "Max rank cutoff (master slider)", 
                min = min(rawdat_res$rawdat$maxRnk, na.rm = TRUE),
                max = max(rawdat_res$rawdat$maxRnk, na.rm = TRUE),
                value = quantile(rawdat_res$rawdat$maxRnk, probs = 20/nrow(rawdat_res$rawdat)))
    
  })  
  
  output$select_sort <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sortcols <- c("pValueLog", "oddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk")
    
    selectInput("select_sort_i", 
                "Select sort column", 
                choices = sortcols,
                selected = "meanRnk")
    
  })  
  
  output$select_userset <- renderUI({
    
    req(rawdat_res$rawdat)
    
    selectInput("select_userset_i", 
                "Filter by user query set", 
                choices = setchoices())
    
  })  
  
  # barplots
    
  # odds ratio
  
  output$slider_oddsratio <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_oddsratio_i",
                "Odds ratio cutoff",
                min = round(min(rawdat_res$rawdat$oddsRatio, na.rm = TRUE), 3),
                max = round(max(inf.omit(rawdat_res$rawdat$oddsRatio), na.rm = TRUE), 3),
                value = round(min(rawdat_res$rawdat$oddsRatio, na.rm = TRUE), 3))
    })
  
  
  output$run_sum <- renderTable({
    

    req(input$select_sort_i)
    
    data.frame(
      x = 
        c("Start time ", 
        "End time ", 
        "Elapsed time ", 
        "Cache ID ", 
        "Regions ",
        "Genome ",
        "Universe ",
        "Database ",
        "LOLAweb commit used "),
      y = 
        c(as.character(rawdat_res$run_sum$start_time),
          as.character(rawdat_res$run_sum$end_time),
          as.character(rawdat_res$run_sum$run_time),
          as.character(rawdat_res$run_sum$cache_name),
          as.character(rawdat_res$run_sum$query_set),
          as.character(rawdat_res$run_sum$genome),
          as.character(rawdat_res$run_sum$universe),
          as.character(rawdat_res$run_sum$region_db),
          as.character(rawdat_res$run_sum$commit)
        )
    , stringsAsFactors = FALSE)

  }, spacing = "s", colnames = FALSE, align = "l")
  
  
  # retrieve job
  shinyqueue::retrieve(con = con, cache_dir = cacheDir, encrypt = TRUE)
  
  scatterplot_input <- function() {
    
    # conditions for handling infinite log pvalues (i.e. pval = 0 due to perfect overlap )
    # also accounting for infinite oddsRatio
    noinfres <- dat()[!is.infinite(dat()$pValueLog) & !is.infinite(dat()$oddsRatio),]
    
    inflogpval <- dat()[is.infinite(dat()$pValueLog),]
    inflogpval$pValueLog <- max(inf.omit(rawdat_res$rawdat$pValueLog))+1
    
    infor <- dat()[is.infinite(dat()$oddsRatio),]
    infor$oddsRatio <- -1e-6
    
    infvals <- rbind(inflogpval, infor)
    
    # case when all the rows have infinite values
    if (all(is.infinite(dat()$pValueLog))) {
      
      p <- 
        ggplot() +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0("Log p-value: ", "Infinite", "\n", "Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ",collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = inflogpval)
      
    } else if (all(is.infinite(dat()$oddsRatio))) {
      p <- 
        ggplot() +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0("Log p-value: ", pValueLog, "\n", "Odds ratio: ", "NA", "\n", "Support: ", support, "\n", "Collection: ",collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = infor)
      
      # ... when some have infinite values
    } else if (any(is.infinite(dat()$pValueLog)) & !any(is.infinite(dat()$oddsRatio))) {
      
      p <-
        ggplot() +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0("Log p-value: ", "Infinite", "\n", "Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ",collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = inflogpval) +
        geom_point(aes(pValueLog, 
                       oddsRatio, 
                       text = paste0("Log p-value: ", pValueLog, "\n", "Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ", collection, "\n", "Description: ", axis_label),
                       col=userSet, 
                       size = log(support)), 
                   alpha=.75, 
                   data = noinfres)
      
      # ... when some have infinite values
    } else if (!any(is.infinite(dat()$pValueLog)) & any(is.infinite(dat()$oddsRatio))) {
      
      p <-
        ggplot() +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0("Log p-value: ", pValueLog, "\n", "Odds ratio: ", "NA", "\n", "Support: ", support, "\n", "Collection: ",collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = infor) +
        geom_point(aes(pValueLog, 
                       oddsRatio, 
                       text = paste0("Log p-value: ", pValueLog, "\n", "Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ", collection, "\n", "Description: ", axis_label),
                       col=userSet, 
                       size = log(support)), 
                   alpha=.75, 
                   data = noinfres) 
      
      # ... when some have infinite values
    } else if (any(is.infinite(dat()$pValueLog)) & any(is.infinite(dat()$oddsRatio))) {
      
      p <-
        ggplot() +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0("Log p-value: ", "Infinite", "\n", "Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ", collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = inflogpval) +
        geom_point(aes(pValueLog, oddsRatio, 
                       # need to construct custom text since the value is fudged
                       text = paste0( "Log p-value: ", pValueLog, "\n", "Odds ratio: ", "NA", "\n", "Support: ", support, "\n", "Collection: ", collection, "\n", "Description: ", axis_label)),
                   col = "black",
                   pch = "O",
                   alpha = 0.75, data = infor) +
        geom_point(aes(pValueLog, 
                       oddsRatio, 
                       text = paste0("Log p-value: ",pValueLog,"\n","Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ", collection,"\n", "Description: ", axis_label),
                       col=userSet, 
                       size = log(support)), 
                   alpha=.75, 
                   data = noinfres)
      
    # ... when none have infinite values
    } else {
      
      p <- 
        ggplot() +
        geom_point(aes(pValueLog, 
                       oddsRatio, 
                       text = paste0("Log p-value: ",pValueLog,"\n","Odds ratio: ", oddsRatio, "\n", "Support: ", support, "\n", "Collection: ", collection,"\n", "Description: ", axis_label),
                       col=userSet, 
                       size = log(support)), 
                   alpha=.75, 
                   data = noinfres)
      
    }
    
    p +         
    xlab("log(p value)") +
    ylab("Odds ratio") +
    scale_size_continuous(range = c(0.5,4)) +
    scale_y_continuous(limits = c(min(rawdat_res$rawdat$oddsRatio)-1, 
                                  max(inf.omit(rawdat_res$rawdat$oddsRatio)))) +
                       # expand = expand_scale(mult = c(0,0))) +
    scale_x_continuous(limits = c(min(rawdat_res$rawdat$pValueLog), 
                                  max(inf.omit(rawdat_res$rawdat$pValueLog))+5)) +
                       # expand = expand_scale(mult = c(0,0))) +
    guides(size = FALSE) +
    theme_ns()
    
  }
  
  output$scatter <- renderPlotly({
    
    req(input$select_collection_i)
    
    plot_render$state <- TRUE
  
    ggplotly(scatterplot_input(), tooltip = "text") %>%
      config(displayModeBar = FALSE) %>%
      layout(
        height = 600,
        width = 600,
        legend = list(orientation = "h",
                      x = 0.44,
                      y = -0.2,
                      title = ""))
    
  })
  
  # download handler
  output$scatterplot_dl <- downloadHandler(
    filename = function() { paste("scatter", 
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = scatterplot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9))
             , device = "pdf")
    }
  )
  
  # call plot
  output$oddsratio_plot <- renderPlot({

    req(input$select_sort_i)

    plot_input(dat(), "oddsRatio", "Odds ratio", input$select_sort_i)

  })
  
  # download handler
  output$oddsratio_plot_dl <- downloadHandler(
    filename = function() { paste("oddsratio", 
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = plot_input(dat(), "oddsRatio", "Odds ratio", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9))
, device = "pdf")
    }
  )
  
  # support
  
  # slider
  output$slider_support <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_support_i",
                "Support cutoff",
                min = round(min(rawdat_res$rawdat$support, na.rm=TRUE), 3),
                max = round(max(rawdat_res$rawdat$support, na.rm=TRUE), 3),
                value = round(min(rawdat_res$rawdat$support, na.rm=TRUE), 3))
    
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
      ggsave(file, plot = plot_input(dat(), "support", "Support", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9))
, device = "pdf")
    }
  )
  
  # pvalue
  
  # slider input
  
  output$slider_pvalue <- renderUI({
    
    req(rawdat_res$rawdat)
    
    sliderInput("slider_pvalue_i", 
                "P-value cutoff", 
                min = round(min(rawdat_res$rawdat$pValueLog, na.rm=TRUE), 3), 
                max = round(max(inf.omit(rawdat_res$rawdat$pValueLog), na.rm=TRUE), 3),
                value = round(min(rawdat_res$rawdat$pValueLog, na.rm=TRUE), 3))
    
    
  })  
  
  output$pvalue_plot <- renderPlot({

    req(input$select_sort_i)

    plot_input(dat(), "pValueLog", "log(p value)", input$select_sort_i)


  })
  
  output$pvalue_plot_dl <- downloadHandler(
    filename = function() { paste("pvalue",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = plot_input(dat(), "pValueLog", "log(p value)", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9))
, device = "pdf")
    }
  )
  
  # genomic distribution plot
  # set up function
  distrib_plot_input <- function() {
    
    genDist <- rawdat_res$genDist
    
    if (is.null(genDist)) {
      
      missing_plot()
      
    } else {
      
      plotChromBins(genomeAggregate = genDist) +
        guides(fill=guide_legend(title="User set"),
               col = guide_legend(title="User set"))
      
    }
    
  }
  
  output$distrib_plot <- renderPlot({

    req(input$select_sort_i)

    distrib_plot_input()

  })
  
  
  # feature distance plot
  dist_plot_input <- function() {
    
    TSSDist <- rawdat_res$TSSDist
    
    if (is.null(TSSDist)) {
      
      missing_plot()
      
    } else {
      
      plotFeatureDist(TSSDist, featureName="TSS") +
        guides(fill=guide_legend(title="User set"),
               col = guide_legend(title="User set"))
      
    }

  }
  
  output$dist_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    dist_plot_input()
    
  })
  
  # partitions plot
  part_plot_input <- function() {
    
    gp <- rawdat_res$gp
    
    if(is.null(gp)) {
      
      missing_plot()
      
    } else {
      
      plotPartitions(gp)
      
    }
  
    
  }
  
  output$part_plot <- renderPlot({
    
    req(input$select_sort_i)
    
    part_plot_input()
    
  })

  # download handler
  output$distrib_plot_dl <- downloadHandler(
    filename = function() { paste("gendist",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = distrib_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)), device = "pdf", width = 11, height = 5)
    }
  )
  
  output$dist_plot_dl <- downloadHandler(
    filename = function() { paste("tssdist",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = dist_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)), device = "pdf")
    }
  )
  
  output$part_plot_dl <- downloadHandler(
    filename = function() { paste("partitions",
                                  ".pdf", 
                                  sep="") },
    content = function(file) {
      ggsave(file, plot = part_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)), device = "pdf")
    }
  )
  
  # to zip all plots and save individually ...
  
  output$all_plots_dl <- downloadHandler(
    filename = function() {
      paste("lolawebplots", "zip", sep=".")
    },
    content = function(fname) {

      ggsave(filename = "scatter.pdf",
             plot = scatterplot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")

      ggsave(filename = "oddsratio.pdf",
             plot = plot_input(dat(), "oddsRatio", "Odds ratio", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")

      ggsave(filename = "support.pdf",
             plot = plot_input(dat(), "support", "Support", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")

      ggsave(filename = "pvalue.pdf",
             plot = plot_input(dat(), "pValueLog", "log(p value)", input$select_sort_i) + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")
      
      # using pdf() device here bc custom aspect ratio not working with ggsave() in this case
      pdf("plots/gendist.pdf", width = 11, height = 5)
      print(distrib_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)))
      dev.off()

      ggsave(filename = "tssdist.pdf",
             plot = dist_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")

      ggsave(filename = "partitions.pdf",
             plot = part_plot_input() + theme(axis.text = element_text(size = 9), text = element_text(size = 9)),
             device = "pdf",
             path = "plots")
      
      # identify files to be zipped and zip them
      fs <- list.files("plots", full.names = TRUE)
      zip(zipfile=fname, files=fs)
      
      # delete tmp plot files after zip is done
      unlink(fs)

    },
    contentType = "application/zip"

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
                                 title = paste0("LOLAweb_", query()[[1]], ".csv"),
                                 extend = "csv",
                                 text = '<i class="fa fa-download"></i> Download CSV')),
                               paging = FALSE)) %>%
      formatRound(columns=c('oddsRatio', 'pValueLog'),
                  digits = 4)
  })
  
}

shinyApp(ui = ui, server = server)