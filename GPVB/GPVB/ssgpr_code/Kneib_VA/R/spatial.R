spatial <- function(geo, map,a_1=1, b_1=.0001){
       if(is.list(map)){gramap <- bnd2gra(map)}else{gramap <- map}
      bsp <-  matrix(0, nrow=length(geo), ncol=ncol(gramap))

for(i in 1:length(geo)){
 bsp[i,which(colnames(gramap)==geo[i])]<-1
}
       attr(bsp, "K")<-gramap
 attr(bsp, "class") <- "geo"
       attr(bsp, "name") <- "geo"
          attr(bsp, "a_1") <-a_1
     attr(bsp, "b_1") <- b_1
attr(bsp, "data")<- geo
     return(bsp)
}