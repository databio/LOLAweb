# misc functions for shinyLOLA app

round_top <- function(x,n) {
  
  round(head(sort(x, decreasing = TRUE), n)[n], 3)
  
}

# infinite value omission

inf.omit <- function(x) {
  
  if(any(is.infinite(x))) {
    
    x <- x[!is.infinite(x)]
    
  }
  
  x
  
}

# plot input

plot_input <- function(res, metric, ylabel, sortcol) {
  
  # limit to 50 results for histograms so they don't get too dense
  if(nrow(res) > 50) {
    
    res <- head(res,50)
  }

  # conditional for inverting rank sorting
  if(grepl("rnk", sortcol, ignore.case = TRUE)) {
    
    # need to order data frame by sort col if it's a rank
    res <- res[order(as.data.frame(res)[,sortcol]), ]
    
    # now construct base layer for plot with reverse on the sort
    p <- ggplot(res, aes(reorder(axis_label, 
      rev(eval(parse(text = sortcol)))), eval(parse(text = metric)), fill = userSet, group = id))
    
  } else {
    
    p <- ggplot(res, aes(reorder(axis_label, 
      eval(parse(text = sortcol))), eval(parse(text = metric)), fill = userSet, group = id))
    
  }

  # Make label size dynamic based on number of entries in the table.
  label_size = 9 + ( 8 * (1 - NROW(res)/50) )
  p +
    geom_bar(stat="identity", position="dodge") +
    xlab("Description") +
    ylab(ylabel) +
    coord_flip() +
    theme_ns() +
    theme(axis.text=element_text(size=label_size)) +
    guides(fill=guide_legend(title="User set"))
  
}

# missing_plot() creates 
missing_plot <- function() {
  
  ggplot(data.frame(x=1:5, y = 1:5), aes(x, y)) +
  geom_blank() +
  annotate("text", x = 3, y = 4, label = "Unable to render plot", size = 10, col = "red") +
  xlab("") +
  ylab("") +
  theme(line = element_blank(), 
        rect = element_blank(), 
        axis.ticks.length = unit(0,"cm"),
        axis.text = element_blank(),
        legend.position = "none", 
        panel.spacing = unit(0,"lines"), 
        plot.margin = unit(c(0, 0, 0, 0), "lines"))
}

# This function just wraps the base Sys.gentenv function to provide a default
# value for the case that the environment variable is not specified.
getEnv = function(envVar, default="") {
  var = Sys.getenv(envVar)
  if (var == "") {
    return(default)
  } else {
    return(var)
  }
}



# Set up cache dir based on environment variable
localDir = getEnv("LWLOCAL", "/data/lola/")
cacheDir = paste0(localDir, "cache/")
logDir = paste0(localDir, "shinylog/")
resultsDir = paste0(localDir, "results/")


refDir = getEnv("LWREF", "/mnt/q/shefflab/LOLAweb/")

dbDir = paste0(refDir, "databases/")
universeDir = paste0(refDir, "universes/")
exampleDir =  paste0(refDir, "examples/")
message("Local dir: ", localDir)
message("Cache dir: ", cacheDir)
message("Reference data dir: ", refDir)
message("universe dir: ", universeDir)
message("dbDir data dir: ", dbDir)
message("exampleDir dir: ", exampleDir)
