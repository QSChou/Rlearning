# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, Ryan J. Tibshirani, and Larry Wasserman. "Nonparametric Modal Regression." arXiv preprint arXiv:1412.1716 (2014).
# Date: 04/07/2015
source("MDR.R")

  ## Generate Data
D1 = ThreeMix()

  ## Mesh points
G.x0 = seq(from=min(D1[,1]), to=max(D1[,1]), length.out = 100)
G.y0 = seq(from = min(D1[,2]), to=max(D1[,2]), length.out=5)
G0 = expand.grid(G.x0,G.y0)


  ## Partial mean shift
h0 = 0.075
D1.pMS = pMS(D1, G0, h0)


par(mar=c(4,4,1,1))
plot(D1,xlab="X",ylab="Y")
points(D1.pMS,cex=0.7,pch=19,col="dodgerblue")


  ## Using "MDR" object
G0.MDR = MDR(D1,G0,h0)

col0 = c("dodgerblue","limegreen")
par(mar=c(4,4,1,1))
plot(D1,xlab="X",ylab="Y", pch=20)
for(j in 1:max(G0.MDR$labels)){
  lines(G0.MDR$ordered[[j]],cex=0.7,pch=19,col=col0[j], lwd=8)
}


  ## Clustering
D1.MDR = MDR(D1, h=h0)

col0 = c("dodgerblue","limegreen", "violet")
par(mar=c(4,4,1,1))
plot(D1,xlab="X",ylab="Y", col=col0[D1.MDR$label], pch=20)


  ## Bandwidth selection
alpha=0.9
h_seq = 1:20/100

CV_error = CV.h(D1, h_seq, alpha=alpha)

plot(CV_error, lwd=5, type="l", ylab="", xlab="h", main=paste("Size of ",100*alpha,"% Prediction interval", sep=""), col="dodgerblue", ylim=c(0,0.5), xlim=c(0,0.2), cex.lab=2, cex.main=2)
mtext("(CV) Size of prediction set", side=2, line=2.2, cex=2)
abline(v= CV_error[which(CV_error[,2]==min(CV_error[,2])),1], col="violet", lwd=4)


  ## Using the optimal bandwidth
h_opt = CV_error[which(CV_error[,2]==min(CV_error[,2])),1]
G0.MDR = MDR(D1,G0,h_opt)
D1.MDR = MDR(D1, h=h_opt)



  ## Prediction set using residual
alpha=0.9
MS_ps = quantile(abs(D1.MDR$res), alpha)

col_ps = c("blue", "dodgerblue")
par(mar=c(4,4,1,1))
plot(D1,xlab="X",ylab="Y", pch=20)
for(j in 1:max(G0.MDR$labels)){
  lines(G0.MDR$ordered[[j]],cex=0.7,pch=19,col=col_ps[1], lwd=8)
  lines(G0.MDR$ordered[[j]]+cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=5)
  lines(G0.MDR$ordered[[j]]-cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=5)
}
legend("bottomleft",c("Modal Regression", paste(100*alpha,"% PS, Modal", sep="")),lwd=c(8,5), col=c("blue","dodgerblue") )



  ## Local regression
span_seq = seq(from=0.1, to=0.9, length.out=100)
loc_PS_seq = rep(NA,100)

for(i in 1:100){
  fit_loc_tmp = loess(D1[,2]~D1[,1],span=span_seq[i])
  loc_PS_seq[i] = quantile(fit_loc_tmp$res, alpha)
}
s_opt = span_seq[which(loc_PS_seq==min(loc_PS_seq))]
#the optimal one

loess_fit  = loess(D1[,2]~D1[,1],span=s_opt)


  ## Prediction set
plot(D1, cex=0.8 , main="Local Regression", pch=20)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)], col="red", lwd=4)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]+quantile(abs(loess_fit$res),alpha), col="orange", lwd=3)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]-quantile(abs(loess_fit$res),alpha), col="orange", lwd=3)
legend("topright",c("Local Regression", paste(alpha,"% PI", sep="")),lwd=c(5,5), col=c("red","orange") )



  ## Comparison of prediction set
plot(D1, cex=0.8 , main=paste("Comparison of ", 100*alpha, "% PS", sep=""), pch=20)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)], col="red", lwd=4)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]+quantile(abs(loess_fit$res),alpha), col="orange", lwd=3)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]-quantile(abs(loess_fit$res),alpha), col="orange", lwd=3)
for(j in 1:max(G0.MDR$labels)){
  lines(G0.MDR$ordered[[j]],cex=0.7,pch=19,col=col_ps[1], lwd=4)
  lines(G0.MDR$ordered[[j]]+cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=3)
  lines(G0.MDR$ordered[[j]]-cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=3)
}
legend("topright",c("Modal Regression", "Local Regression", paste(100*alpha,"% PS, Modal", sep=""),paste(100*alpha,"% PS, Local", sep="")),lwd=c(5,5,5,5), col=c("blue","red","dodgerblue","orange") )


  ## Size of prediction set
plot(CV_error, lwd=5, type="l", ylab="", xlab="h", main="", cex.lab=2, col="dodgerblue")
abline(h= quantile(abs(loess_fit$res),alpha), col="red", lwd=5)
mtext("(CV) Size of prediction set", side=2, line=2.2, cex=2)
legend("topright", c("Modal Regression","Local Regression"), col=c("black","red"), lwd=c(5,5))



  ## Mixture regression
library(mixtools)

mix1 = regmixEM.lambda(D1[,2],D1[,1], k=3)
for(i in 1:10000){
  mix_tmp = regmixEM.lambda(D1[,2],D1[,1], k=3)
  if(mix1$loglik>mix_tmp$loglik)
  {
    mix1 = mix_tmp	
  }
}	##find the optimalEM over several runs


  ## Comparison of prediction sets
plot(D1, cex=0.8, pch=20)
for(i_k in 1:3){
  abline(a=mix1$beta[1,i_k], b=mix1$beta[2,i_k], lwd=4, col="purple")
  abline(a=mix1$beta[1,i_k]+1.96*mix1$sigma[i_k], b=mix1$beta[2,i_k], lwd=3, col="violet")
  abline(a=mix1$beta[1,i_k]-1.96*mix1$sigma[i_k], b=mix1$beta[2,i_k], lwd=3, col="violet")
}
for(j in 1:max(G0.MDR$labels)){
  lines(G0.MDR$ordered[[j]],cex=0.7,pch=19,col=col_ps[1], lwd=4)
  lines(G0.MDR$ordered[[j]]+cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=3)
  lines(G0.MDR$ordered[[j]]-cbind(rep(0, nrow(G0.MDR$ordered[[j]])),MS_ps),cex=0.7,pch=19,col=col_ps[2], lwd=3)
}
legend("topright",c("Modal Regression", "Mixture Regression", paste(100*alpha,"% PS, Modal", sep=""),paste(100*alpha,"% PS, Mixture", sep="")),lwd=c(5,5,5,5), col=c("blue","purple","dodgerblue","violet") ,bg="white")

