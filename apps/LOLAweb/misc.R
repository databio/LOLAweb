# misc functions for shinyLOLA app

round_top <- function(x,n) {
  
  round(head(sort(x, decreasing = TRUE), n)[n], 3)
  
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
    dat <- res[order(as.data.frame(res)[,sortcol]), ]
    
    # now construt base layer for plot with reverse on the sort
    p <- ggplot(dat, aes(reorder(axis_label, rev(eval(parse(text = sortcol)))), eval(parse(text = metric)), fill = userSet, group = id))
    
  } else {
    
    p <- ggplot(res, aes(reorder(axis_label, eval(parse(text = sortcol))), eval(parse(text = metric)), fill = userSet, group = id))
    
  }
  
  p +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Description") +
    ylab(ylabel) +
    coord_flip() +
    theme_ns()
  
}