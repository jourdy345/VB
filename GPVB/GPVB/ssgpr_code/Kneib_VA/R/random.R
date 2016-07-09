random <- function(ind, a_1=1, b_1=1e-5){
  m <- length(table(ind))
  n <- length(ind)
  V <- table(which(ind==ind),ind)
  attr(V, "K")<- diag(1, m)
  attr(V,"a_1")<-a_1
  attr(V,"b_1")<-b_1
  attr(V, "class") <- "random"
  attr(V, "name") <- "random" 
  return(V)
  
}