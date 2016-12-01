library(R.matlab)

#setwd("C:/Users/Physium/Dropbox/Edge Selection/Work with David/VBILL code/Simple Normal/normal results")
setwd("~/Dropbox/Edge Selection/Work with David/VBILL code/Simple Normal/normal results")
source("function.R")

#n is sample size of data
#N is number of summary statistics for the VBIl/SBL approach

SBL4 <- readMat("SBL4var_50_Anal_QMC.mat")  # SBL N = 50 using analytic, n = 4
SBL4I <- readMat("SBL4var_50_Monte_QMC.mat") # SBL N = 50 using monte carlo, n = 4
SBL8I <- readMat("SBL8var_50_Monte_QMC.mat") # SBL N = 50 using monte carlo, n = 8

SBL8I_normal <- readMat("SBL8var_50_Monte.mat")

SBL4I_100 <- readMat("SBL4var_100_Monte_QMC.mat") # SBL n = 100 using monte carlo, n = 4
SBL8I_100 <- readMat("SBL8var_100_Monte_QMC.mat") # SBL n = 100 using monte carlo, n = 8

VBIL4smallI <- readMat("VBIL4var01_Monte_QMC.mat") # VBIL, V(log likelihood) < 0.1, using monte carlo, n = 4
VBIL4large <- readMat("VBIL4var05_Anal_QMC.mat") # VBIL, V(log likelihood) < 0.5, using analytic, n = 4
VBIL4largeI <- readMat("VBIL4var05_Monte_QMC.mat") # VBIL, V(log likelihood) < 0.5, using monte carlo, n = 4

VBIL8large_normal <- readMat("VBIL8var05_Anal.mat")
VBIL8largeI_normal <- readMat("VBIL8var05_Monte.mat")

VBIL8smallI <- readMat("VBIL8var01_Monte_QMC.mat") # VBIL, V(log likelihood) < 0.1, using monte carlo, n = 8
VBIL8largeI <- readMat("VBIL8var05_Monte_QMC.mat") # VBIL, V(log likelihood) < 0.5, using monte carlo, n = 8


####### This part tries to compare the difference between analytic and monte carlo estimates for VBSL ######

pdf("SBL4.pdf",width=9,height=5)

mtest <- matrix(c(1,1,1,2,3,4),nrow = 2,ncol = 3,byrow = TRUE)

layout(mat = mtest,heights = c(1,10))
  par(mar = c(0,0,0,0))
plot_colors <- c(1,4,2)
text <- c("True posterior", "VBSL (Analytic)","VBSL")
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.5*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(3,1,1), lwd=2, cex=1.5, horiz = TRUE,seg.len=4)
par(xpd=FALSE)
  par(mar =c (7, 5, 3, 2))

Normal_post(SBL4,SBL4I,title="")
Normal_grad(SBL4,SBL4I,title="")
Normal_LB(SBL4,SBL4I,title="",1)

dev.off()

####### This part tries to compare the difference between analytic and monte carlo estimates for VBIL ######

pdf("VBIL4.pdf",width=9,height=5)

mtest <- matrix(c(1,1,1,2,3,4),nrow = 2,ncol = 3,byrow = TRUE)
layout(mat = mtest,heights = c(1,10))
  par(mar = c(0,0,0,0))
plot_colors <- c(1,4,2)
text <- c("True posterior", "VBIL (Analytic)","VBIL")
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.5*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(3,1,1), lwd=2, cex=1.5, horiz = TRUE,seg.len=4)
par(xpd=FALSE)
  par(mar =c (7, 5, 3, 2))

Normal_post(VBIL4large,VBIL4largeI,ylim=c(0,1),title="")
Normal_grad(VBIL4large,VBIL4largeI,title="")
Normal_LB_NULL(VBIL4large,VBIL4largeI,title="",2)
dev.off()

######### VBSL and VBIL under different setting of N for n = 4 ##########

pdf("Normal_4_N.pdf",width=9,height=5)
mtest <- matrix(c(1,1,2,2,3,4,5,6),nrow = 2,ncol = 4,byrow = TRUE)

layout(mat = mtest,heights = c(1,5))
  par(mar = c(0,0,0,0))
plot_colors <- c(4,2)
text <- c("VBSL (N=50)","VBSL (N=100)")
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.2*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(1), lwd=2, cex=1.2, horiz = TRUE,seg.len=3)
par(xpd=FALSE)

#layout(mat = mtest,heights = c(1,5))
  par(mar = c(0,0,0,0))
plot_colors <- c(4,2)
text <- c(expression(paste(VBIL[0.1])),expression(paste(VBIL[0.5])))
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.2*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(1), lwd=2, cex=1.2, horiz = TRUE,seg.len=3)
par(xpd=FALSE)

  par(mar =c (5, 4, 2, 2))

Normal_post(SBL4I,SBL4I_100,title="")
Normal_LB(SBL4I,SBL4I_100,title="",1)


Normal_post(VBIL4smallI,VBIL4largeI,ylim=c(0,1),title="")
Normal_LB(VBIL4smallI,VBIL4largeI,title="",2)
dev.off()

######### VBSL and VBIL under different setting of N for n = 8 ##########

pdf("Normal_8_N.pdf",width=9,height=5)
mtest <- matrix(c(1,1,2,2,3,4,5,6),nrow = 2,ncol = 4,byrow = TRUE)

layout(mat = mtest,heights = c(1,5))
  par(mar = c(0,0,0,0))
plot_colors <- c(4,2)
text <- c("VBSL (N=50)","VBSL (N=100)")
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.2*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(1), lwd=2, cex=1.2, horiz = TRUE,seg.len=3)
par(xpd=FALSE)

#layout(mat = mtest,heights = c(1,5))
  par(mar = c(0,0,0,0))
plot_colors <- c(4,2)
text <- c(expression(paste(VBIL[0.1])),expression(paste(VBIL[0.5])))
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = 1.2*max(sapply(text, strwidth)),
       col=plot_colors, lty=c(1), lwd=2, cex=1.2, horiz = TRUE,seg.len=3)
par(xpd=FALSE)

  par(mar =c (5, 4, 2, 2))

Normal_post(SBL8I,SBL8I_100,title="")
Normal_LB(SBL8I,SBL8I_100,title="",1)


Normal_post(VBIL8smallI,VBIL8largeI,ylim=c(0,1.7),title="")
Normal_LB(VBIL8smallI,VBIL8largeI,title="",2)
dev.off()


