coef.va <- function(object, ...){
  x <- object
  co <- as.vector(x$be)
names(co)<- rownames(x$be)
  return(co)
}