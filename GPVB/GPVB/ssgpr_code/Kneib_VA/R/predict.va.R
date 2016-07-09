predict.va<- function(object, newdata=NULL, ...){
   L <- object
   if(is.null(newdata)){
    prediction <-L$prediction
      }else{
    if(length(L$gam)>0){spli=T}else{spli=F}
    if(spli){
       if("random"%in%names(L$gamma)){
         
         L$gamma[[which(names(L$gamma)=="random")]]<-rep(0, length(L$gamma[[which(names(L$gamma)=="random")]]))}
       #print(L$gamma)
            #mu_gamG <-as.vector(lapply(L$gam, rbind))
    }

    
    if (length(L$bet)>0){li <- TRUE
                          mu_bet <- L$bet
                          if(rownames(L$bet)[1]=="intercept"){
                            intercept=TRUE
                          }else{intercept=FALSE}}else{li <- FALSE}
  formula <- L$call
  datapred <- as.data.frame(newdata)  
  Zpred <- list()
  j <- 1
  n <- nrow(datapred)
Xlinearpred <- as.matrix(rep(1,n))
if(length(which(colnames(datapred)==as.character((formula)[2])))==0){
  datapred<- cbind(1, datapred)
  colnames(datapred)[1]<-as.character((formula)[2])
}

type <- vector(length = length(model.frame(formula, datapred))-1)
for(l in 2:length(model.frame(formula, data=datapred))){
  b <- model.frame(formula, datapred)[,l]
    if(is.null(attr(b,"class"))){Xlinearpred<-cbind(Xlinearpred,b)}  
  else{if(attr(b,"class")%in%c("bspline", "geo", "random")){
    if(attr(b,"class")=="bspline"){
    pas<-which(names(L$Z)== names(get_all_vars(formula)[l]))
    knots<-L$knots[[pas]]
    mx <- L$mx[[pas]]
    b <- bsplines(attr(b, "data"), knots=knots, mx=mx)}
    
    Zpred[[j]]<-as.matrix(b)
  j<-j+1}
  }
    
  }
  

  q=length(Zpred)
  
  if(spli){

   Zpredg <- Zpred[[1]]
   mu_gamG<-L$gam[[1]]
   if(q >1){
   for(j in 2:q){
 	     Zpredg <- cbind(Zpredg, Zpred[[j]])
      mu_gamG<-c(mu_gamG,L$gam[[j]])
   }
   }
}
if(intercept==FALSE){Xlinearpred <- Xlinearpred[,-1]}
if(spli&li) {
 
prediction <- Zpredg%*%mu_gamG + Xlinearpred%*%mu_bet}else{
  if(spli){
    prediction <- Zpredg%*%mu_gamG     
  }else{
    if(li){
      prediction <-  Xlinearpred%*%mu_bet
    }
  }
}
}
  return(prediction)}