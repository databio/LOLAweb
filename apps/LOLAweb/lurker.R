# get environment variables captured in misc.R
source("misc.R")
library(LOLA)
library(GenomicRanges)
library(simpleCache)
library(sodium)
library(GenomicDistributions)

lwcrunch <- function() {
  
  # get lw version
  lw_version <- system(command = "git rev-parse HEAD | cut -c1-9", 
                       intern = TRUE)
  
  refgenome <- unlist(input$refgenome)
  
  message("Calculating region set enrichments ...")
  
  start_time <- Sys.time()
  
  run_time <- system.time({
    
    # user set(s)
    userSets <- list()
    
    # if(!exampleuserset$toggle) {
    # if(!is.null(input$userset)) {
    
    input_userset  <- as.data.frame(input$userset)
    
    if(ncol(input_userset) != 0) {
        
      message("Loading uploaded data")
    
      
      # for (i in 1:length(input$userset[,1])) {
      for (i in 1:nrow(input_userset)) {
          
        
        # fixing file path to env var for tempdir
        fp <- paste0(tempDir, input_userset[i, 'name'])
        
        # userSet <- readBed(input_userset[i, 'datapath'])
        userSet <- readBed(fp)
        
        userSets[[i]] <- userSet
        
      }
      
      #print(list(userSet,userSets))
    
      userSets = GRangesList(userSets)
      
      names(userSets) = input_userset[,"name"]
      
    } else {
      
      message("Loading example data")
      
      datapath <- paste0(exampleDir, refgenome, "/", input$defaultuserset)
      
      userSet = readBed(datapath)
      userSets[[1]] <- userSet
      
      userSets = GRangesList(userSets)
      
      names(userSets) = input$defaultuserset
      
    }
    
    # universe
    
    if(input$universe_opts == "user") {
      
      userUniverse <- readBed(file = input$useruniverse$datapath)
      
      universename <- input$useruniverse$datapath
      
    } else if (input$universe_opts == "default") {
      
      datapath <- paste0(universeDir, 
                         refgenome, "/", 
                         unlist(input$defaultuniverse))
      
      userUniverse <- readBed(file = datapath)
      
      universename <- input$defaultuniverse
      
    } else if (input$universe_opts == "build") {
      
      userUniverse <- buildRestrictedUniverse(userSets)
      
      universename <- paste0(c(names(userSets), 
                               "universe"),
                             collapse = "_")
      
    }
    
    userSetsRedefined =	redefineUserSets(userSets, userUniverse)
    
    # set up path to databases for given reference genome
    dbPath = paste(dbDir,
                   input$loladb,
                   refgenome,
                   sep = "/"
    )
    
    # load region data for the given reference genome
    regionDB = loadRegionDB(dbLocation=dbPath)
    
    # garbage collection
    gc()
    
    resRedefined = runLOLA(userSetsRedefined,
                           userUniverse,
                           regionDB,
                           cores=4)
    
    # need to make sure user set is discrete even if coded as number
    resRedefined$userSet = as.character(resRedefined$userSet)
    
    gd_data <- utils::data(package="GenomicDistributions")$results[,"Item"]
    
    # calculate distribution over chromosomes for plotting
    if (paste0("chromSizes_", refgenome) %in% gd_data) {
      
      genDist = aggregateOverGenomeBins(userSets, refgenome)
      
    } else {
      
      genDist = NULL
      
    }
    
    # calculate distances to TSSs
    if (paste0("TSS_", refgenome) %in% gd_data) {
      
      TSSDist = TSSDistance(userSets, refgenome)
      
      
    } else {
      
      TSSDist = NULL
      
    }
    
    # distribution of overlaps for a query set to genomic partitions
    if (paste0("geneModels_", refgenome) %in% gd_data) {
      
      gp = genomicPartitions(userSets, refgenome)
      
    } else {
      
      gp = NULL
      
    }
    
  })
  
  run_sum <- 
    data.frame(
      start_time = as.character(start_time),
      end_time = as.character(Sys.time()),
      run_time = paste0(round(run_time[3])+1, " seconds"),
      cache_name = unlist(input$job_id),
      query_set = paste(gsub(".bed","",unique(resRedefined$userSet)), collapse = "\n"),
      genome = refgenome,
      universe = gsub(".bed", "", universename),
      region_db = unlist(input$loladb),
      commit = lw_version
    )
  
  # create named list of multiple objects for plotting
  res = list(resRedefined = resRedefined,
             genDist = genDist,
             TSSDist = TSSDist,
             run_sum = run_sum,
             gp = gp)
  
  # caching
  # need to call keyphrase from reactive above because it is used to construct link
  key <- hash(charToRaw(unlist(input$job_id)))
  msg <- serialize(res,connection = NULL)
  
  cipher <- data_encrypt(msg, key)
  
  simpleCache(cacheName = unlist(input$job_id), 
              instruction = { cipher },
              noload = TRUE)
  
  return(res)
  
}

# what is the name of the db for the jobs
db_url <- "mongodb://localhost:27017"
db_name <- "lolaweb"

# # where do the caches live?
# cache_dir <- "cache/"

# establish connection to jobdb
con <- shinyqueue::connect(db_url = db_url, db_name = db_name)

shinyqueue::lurk(process = list("lolaweb" = lwcrunch()),
                                con = con)
