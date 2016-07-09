risk <-
function(y, pred,quant=FALSE, tau=0.5){
  u <- y-pred
  if(quant){
risk<- mean(u * (tau - (u<0)))}else{
  risk<- mean(u^2)
  }
return(risk)    
}

