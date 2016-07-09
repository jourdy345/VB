lin<- function(x, mu0=0, sig0=10000){
  

 lin <- as.matrix(x)
  attr(lin, "class") <- "linear"
     attr(lin, "mu0")<- mu0
     attr(lin, "sig0") <-sig0
          return(lin)
}