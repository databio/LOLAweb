# misc functions for shinyLOLA app

round_top <- function(x,n) {
  
  round(head(sort(x, decreasing = TRUE), n)[n], 3)
  
}

# create user set choices to be either multiple or single choice
setchoices <- function() {
  
  req(input$run)
  
  if(length(unique(raw_dat()$userSet)) == 1) {
    
    unique(raw_dat()$userSet)
    
  } else {
    
    c("All User Sets", unique(raw_dat()$userSet))
    
  }
}