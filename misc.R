# misc functions for shinyLOLA app

round_top <- function(x,n) {
  
  round(head(sort(x, decreasing = TRUE), n)[n], 3)
  
}