colfunc <- function(i,gd){
  index = c(i-1,i+1,i-gd,i+gd)
  return(index)
}