summary.va <- function(object, ...){
  x <- object
  V <- list()
  spatial=FALSE
  nonlin=FALSE
  lin=FALSE
  random=FALSE
  V$iterations <- x$i
  V$call <- x$call
if(length(x$tau)){  V$tau <- x$tau}
if(length(x$be)>0){lin=TRUE}
 if(lin){
  M <- matrix(nrow=length(x$be), ncol=4)
  M[,1]<- x$be
  M[,2]<-sqrt(diag(x$sig_be))
  M[,3]<-qnorm(.025,x$be,sqrt(diag(x$sig_be)))
  M[,4]<-qnorm(.975,x$be,sqrt(diag(x$sig_be)))
  colnames(M)<- c("Estimates", "Std. Error", "2.5%-Quantile","97.5%-Quantile")
  rownames(M)<-rownames(x$be)
  V$coefficients <- M
   }
 
 if(length(which(names(x$gam)%in%c("geo", "random")==FALSE))>0){ 
   ind<- which(names(x$gam)%in%c("geo", "random")==FALSE)
   l <- length(ind)
   V$nl <- list()
  for(i in 1:l){
   V$nl[[i]] <- (x$gam)[[ind[i]]]
   names(V$nl)[i]<-names(x$gam)[ind[i]]}
 }
 
 if(length(which(names(x$gam)=="geo" ))>0){
 l <- length(which(names(x$gam)=="geo" )>0)    
 
        qua <- list()
   qua[[1]]<-(1:length(x$gam[[l]]))
   qua[[2]] <- rownames(x$gramap)
   qua[[3]] <-as.matrix(x$gam[[l]])
   
V$geo <- (qua)
 attr(V$geo, "number")<- length(table(x$geo))

   }
   
   if(length(which(names(x$gam)=="random"))>0){
    V$random <- x$gam[[which(names(x$gam)=="random")]] 
    }  
   class(V)<- "summary.va"
  return(V)
 }