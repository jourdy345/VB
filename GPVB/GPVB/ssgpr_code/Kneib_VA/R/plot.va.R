plot.va <- 
  function(x, map=NULL, ask=TRUE, linear = TRUE, bsplines = TRUE, spatial = TRUE, 
           random = TRUE, ...){
   if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
L <- x  
### LINEAR
   if(linear==TRUE){
  if(length(L$be)>0){
  l<- length(L$be)
  if(l<4){
  par(mfrow=c(1,l))}else{
   if(l<10){
   par(mfrow=c(round(sqrt(l)), ceiling(sqrt(l))))}else{par(mfrow=c(3,3))}
  }
  for(i in 1:l){
   m <- L$be[i]
   v <- sqrt(L$sig_bet[i,i])
   span <- seq(qnorm(.0001,m,(v)),qnorm(.9999,m,v), by=v/100)#seq(m-abs(v/m),m+abs(v/m),abs(m*v/100))


   plot(span,dnorm(span,m,(v)), type="l",
   main = paste("density of linear effect of", rownames(L$be)[i]),
   xlab="", ylab="",...)
   q25 <- qnorm(.025, m, v)
   q975<-qnorm(.975, m, v)
   segments(q25, 0, q25, dnorm(q25, m,v) )
   if(0<q975&q25<0){
   segments(0, 0, 0, dnorm(0, m,v), col="red" )}
   segments(q975, 0, q975, dnorm(q975, m,v) )

  }
  }
   }

### B-SPLINES  
  if(bsplines==TRUE){
   if(length(which(names(L$gam)%in%c("geo", "random")==FALSE))){
  ind<- which(names(L$gam)%in%c("geo", "random")==FALSE)
   l <- length(ind)
 if(l<4){
  par(mfrow=c(1,l))}else{
       if(l<10){
   par(mfrow=c(round(sqrt(l)), ceiling(sqrt(l))))}else{par(mfrow=c(3,3))}
  }
  
     
    for(i in ind){
 x <- as.vector(as.matrix(L$Xnonlinear)[,i])
   s <- order(x)
 q25<- qnorm(.025,0,sqrt(diag(L$sig_gamma[[i]])))
 c1 <-(L$Z[[i]]%*%(L$gam[[i]]+q25))[s]
q975<- qnorm(.975,0,sqrt(diag(L$sig_gamma[[i]])))
 c2 <-(L$Z[[i]]%*%(L$gam[[i]]+q975))[s]
 plot(x[s], (L$Z[[i]]%*%L$gam[[i]])[s], 
        main=paste("Nonlinear effect of ",names(L$gam)[i] ),
       ylab = " ", xlab=" ",
       type="l", ylim = c(min(c1,c2), max(c1,c2)),...)
 
 lines(x[s], c1, lty=2)

 lines(x[s],c2 , lty=2)
}
    }
 
}
### GEO
   if(spatial ==TRUE){
 if(length(which(names(L$gam)=="geo" ))>0){
   par(mfrow=c(1,3))
   l <- which(names(L$gam)=="geo" )
  qua <- list()#cbind(1:length(L$gam[[l]]), rownames(L$gramap), as.matrix(L$gam[[l]]))
   qua[[1]]<-(1:length(L$gam[[l]]))
   qua[[2]] <- rownames(L$gramap)
   qua[[3]] <-as.matrix(L$gam[[l]])
   drawmap(as.data.frame(qua), map=map, main="Spatial effect", ...)
   qua2<-list()
   qua2[[1]]<-(1:length(L$gam[[l]]))
   qua2[[2]] <- rownames(L$gramap)
   uq25 <- as.matrix(L$gam[[l]]+qnorm(.025,0,sqrt(diag(L$sig_gamma[[l]]))))
   oq975 <- as.matrix(L$gam[[l]]+qnorm(.975,0,sqrt(diag(L$sig_gamma[[l]]))))
   qua2[[3]] <- rep(0, length(L$gam[[l]]))
   qua2[[3]][which(uq25>0)]<-1
   qua2[[3]][which(oq975<0)]<--1
   drawmap(as.data.frame(qua2), map=map, main="95%-Significance", pcat=T, ...)
   qua3<-list()
   qua3[[1]]<-(1:length(L$gam[[l]]))
   qua3[[2]] <- rownames(L$gramap)
   uq1 <- as.matrix(L$gam[[l]]+qnorm(.1,0,sqrt(diag(L$sig_gamma[[l]]))))
   oq9 <- as.matrix(L$gam[[l]]+qnorm(.9,0,sqrt(diag(L$sig_gamma[[l]]))))
   qua3[[3]] <- rep(0, length(L$gam[[l]]))
   qua3[[3]][which(uq1>0)]<-1
   qua3[[3]][which(oq9<0)]<--1
   drawmap(as.data.frame(qua3), map=map, main="80%-Significance", pcat=T,...)
   }
   }

### RANDOM
   if(random==TRUE){
  if(length(which(names(L$gam)=="random"))>0){
    par(mfrow=c(1,1))
    
    barplot(L$gam$random, space=0.2,...)
    if(length(L$gam$random)<31){
    vek<-1:length(L$gam$random)+seq(.2,.2*length(L$gam$random), length=length(L$gam$random))
    pos <- which(L$gam$random>0)
    neg <- which(L$gam$random<0)
    xnames= names(L$Z$random)
    text(x=vek[pos]-.7, quantile(L$gam$random[neg],.9), cex=1,  xpd=TRUE,xnames[pos], srt=45) 
    text(x=vek[neg]-.3, quantile(L$gam$random[pos], .1), cex=1,  xpd=TRUE,xnames[neg], srt=45)}
}
}
  par(mfrow=c(1,1))
  on.exit(devAskNewPage(FALSE))
}

