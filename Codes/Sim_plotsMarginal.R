## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
setwd('D:/Study/My projects/Depth-model-selection/Codes')
#source('misc_functions.R')

load('lmOutMarginal.rda')

OutList = lmOutMarginal
IC.mat = cbind(OutList[[1]]$IC.vec[,1], OutList[[2]]$IC.vec[,1],
               OutList[[3]]$IC.vec[,1], OutList[[4]]$IC.vec[,1])

par(mfrow=c(2,2))
for(i in 1:4){
  #pdf(paste0('simplot',2*i,'.pdf'),4,4)
  plot(OutList[[i]]$Cn.mat[,1], OutList[[i]]$Cn.mat[,2], type='l', lwd=2,
       xlab="Bootstrap standard deviation", ylab="P(correct model)",
       main=paste0("k=",2*i), col="blue", ylim=c(0,1))
  #lines((1:20)/60, Corr.mat[,i], type='b', col="blue", lwd=2)
  abline(h=IC.mat[1,i], col="red", lty=2, lwd=2)
  abline(h=IC.mat[2,i], col="red", lwd=2)
  abline(h=IC.mat[3,i], lty=2, lwd=2)
  abline(h=IC.mat[4,i], lwd=2)
  #dev.off()
}
par(mfrow=c(1,1))

par(mfrow=c(2,2))
for(i in 1:4){
  #pdf(paste0('simplot',2*i,'bias.pdf'),5,5)
  plot(OutList[[1]]$Cn.mat[,1], OutList[[i]]$Cn.mat[,4], type='l', lwd=2,
       xlab="Bootstrap standard deviation", ylab="absolute bias",
       main=paste0("k=",2*i), col="blue", ylim=c(0,0.1))
  for(j in 5:13){
    if(j-3 <= 2*i){
      lines(OutList[[1]]$Cn.mat[,1], OutList[[i]]$Cn.mat[,j], lwd=2, col="blue")
    }
    else{
      lines(OutList[[1]]$Cn.mat[,1], OutList[[i]]$Cn.mat[,j], lwd=2, lty=2)
    }
  }
  #dev.off()
}
par(mfrow=c(1,1))

pdf('simplot24.pdf',7,4.5)
defaultPar = par()
par(mfrow=c(1,2), oma=c(3,0,1,0))

plot(OutList[[1]]$Cn.mat[,1], OutList[[1]]$Cn.mat[,2], type='l', lwd=2,
     xlab="Bootstrap standard deviation", ylab="P(correct model)",
     main=paste0("k=2"), col="blue", ylim=c(0,1))
abline(h=IC.mat[1,1], col="red", lty=2, lwd=2)
abline(h=IC.mat[2,1], col="red", lwd=2)
abline(h=IC.mat[3,1], lty=2, lwd=2)
abline(h=IC.mat[4,1], lwd=2)

plot(OutList[[1]]$Cn.mat[,1], OutList[[2]]$Cn.mat[,2], type='l', lwd=2,
     xlab="Bootstrap standard deviation", ylab="P(correct model)",
     main=paste0("k=4"), col="blue", ylim=c(0,1))
abline(h=IC.mat[1,2], col="red", lty=2, lwd=2)
abline(h=IC.mat[2,2], col="red", lwd=2)
abline(h=IC.mat[3,2], lty=2, lwd=2)
abline(h=IC.mat[4,2], lwd=2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',
       c("AIC backward","AIC all subset","BIC backward","BIC all subset","Depth-based"),
       col=c("red", "red", "black", "black","blue"),
       lty=c(2,1,2,1,1),
       lwd=2, xpd=T, ncol=3, bty="n")

par(defaultPar)
dev.off()

pdf('simplot68.pdf',7,4.5)
defaultPar = par()
par(mfrow=c(1,2), oma=c(3,0,1,0))

plot(OutList[[1]]$Cn.mat[,1], OutList[[3]]$Cn.mat[,2], type='l', lwd=2,
     xlab="Bootstrap standard deviation", ylab="P(correct model)",
     main=paste0("k=6"), col="blue", ylim=c(0,1))
abline(h=IC.mat[1,3], col="red", lty=2, lwd=2)
abline(h=IC.mat[2,3], col="red", lwd=2)
abline(h=IC.mat[3,3], lty=2, lwd=2)
abline(h=IC.mat[4,3], lwd=2)

plot(OutList[[1]]$Cn.mat[,1], OutList[[4]]$Cn.mat[,2], type='l', lwd=2,
     xlab="Bootstrap standard deviation", ylab="P(correct model)",
     main=paste0("k=8"), col="blue", ylim=c(0,1))
abline(h=IC.mat[1,4], col="red", lty=2, lwd=2)
abline(h=IC.mat[2,4], col="red", lwd=2)
abline(h=IC.mat[3,4], lty=2, lwd=2)
abline(h=IC.mat[4,4], lwd=2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',
       c("AIC backward","AIC all subset","BIC backward","BIC all subset","Depth-based"),
       col=c("red", "red", "black", "black","blue"),
       lty=c(2,1,2,1,1),
       lwd=2, xpd=T, ncol=3, bty="n")

par(defaultPar)
dev.off()
