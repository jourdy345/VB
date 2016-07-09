bsplines<- 
  function(x, M = 20, degree = 3, order = 2, a_1=1, b_1=1e-5, knots=NULL, mx=NULL){
name<- names(as.data.frame(x))

xm <- x
if(is.null(mx)){
mx <- mean(x)}
    x <- x-mx
if(is.null(knots)){
    nrknots <- M
    minx <- min(x)-.001
    maxx <- max(x)+.001
    step <- (maxx-minx)/(nrknots-1)
    inner.knots <- seq(minx,maxx,length=nrknots)
    knots <- seq(minx-3*step,maxx+3*step,by=step)}

 bsp<- spline.des(knots, x, ord=4, outer.ok=TRUE)$design
    nbasis <- M+degree-1
      D <- diag(nbasis)
      for(i in 1:order){
         D <- diff(D)
      tD <- t(D)
      }
       K  <- tD %*% D
attr(bsp, "name") <- name
  attr(bsp, "class") <- "bspline"
     attr(bsp, "K")<- K
     attr(bsp, "a_1") <-a_1
     attr(bsp, "b_1") <- b_1
    attr(bsp, "data")<- xm
attr(bsp, "knots")<-knots
attr(bsp, "mx")<-mx
     return(bsp)

}

