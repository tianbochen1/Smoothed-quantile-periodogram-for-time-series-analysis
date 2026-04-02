library(fields)
library(quantreg)
library(expectreg)
source('fn.r')
library(fGarch)
#######################plot loss
load('GT.rdata')
KK=1000

par(mar=c(3,3,1.5,1.5),mgp=c(3,0.5,0)) 
set.panel(2,2)
image((1:399)/800, y = tau, z=GT[[4]],xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), 
      main=expression('ARMA'))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(Freq),ylab='Quantile', line=2, cex.lab=1.2)

image((1:399)/800, y = tau, z=GT[[8]],xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), 
      main=expression('QAR'))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(Freq),ylab='Quantile', line=2, cex.lab=1.2)

image((1:399)/800, y = tau, z=GT[[12]],xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), 
      main=expression('ARFIMA'))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(Freq),ylab='Quantile', line=2, cex.lab=1.2)

image((1:399)/800, y = tau, z=GT[[16]],xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), 
      main=expression('GARCH'))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab=expression(Freq),ylab='Quantile', line=2, cex.lab=1.2)


########################### 
set.seed(1)
load('GT.rdata')
KK = 1000

dist_arma100 = 0;dist_arma200 = 0;dist_arma400 = 0;dist_arma800 = 0;
dist2_arma100 = 0;dist2_arma200 = 0;dist2_arma400 = 0;dist2_arma800 = 0;


for(i in 1:KK){  
  data1 = as.numeric(arima.sim(list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488),sd=1),n=100))
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 = qpsmo(data1,freq,tau, 1)
  
  data2 = as.numeric(arima.sim(list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488),sd=1),n=200))
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q2 = qpsmo(data2,freq,tau, 1)
  
  data3 = as.numeric(arima.sim(list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488),sd=1),n=400))
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q3 = qpsmo(data3,freq,tau, 1)
  
  data4 = as.numeric(arima.sim(list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488),sd=1),n=800))
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q4 = qpsmo(data4,freq,tau, 1)
  
  for(j in 1:46){
    q1[,j] = smooth.spline(q1[,j])$y
    q2[,j] = smooth.spline(q2[,j])$y
    q3[,j] = smooth.spline(q3[,j])$y
    q4[,j] = smooth.spline(q4[,j])$y
  }
  
  q1 = spec.normalize(q1)
  q2 = spec.normalize(q2)
  q3 = spec.normalize(q3)
  q4 = spec.normalize(q4)
  
  q1[q1 < 0] <- 0
  q2[q2 < 0] <- 0
  q3[q3 < 0] <- 0
  q4[q4 < 0] <- 0
  
  dist_arma100 = dist_arma100 + JSD(q1,GT[[1]])/KK
  dist_arma200 = dist_arma200 + JSD(q2,GT[[2]])/KK
  dist_arma400 = dist_arma400 + JSD(q3,GT[[3]])/KK
  dist_arma800 = dist_arma800 + JSD(q4,GT[[4]])/KK
  
  dist2_arma100 = dist2_arma100 + rmse(q1,GT[[1]])/KK
  dist2_arma200 = dist2_arma200 + rmse(q2,GT[[2]])/KK
  dist2_arma400 = dist2_arma400 + rmse(q3,GT[[3]])/KK
  dist2_arma800 = dist2_arma800 + rmse(q4,GT[[4]])/KK
  if(i %% 20 == 0){
    print(c('arma',i))}
}
####################################  QAR ###  #######################################
dist_qar100 = 0;dist_qar200 = 0;dist_qar400 = 0;dist_qar800 = 0;
dist2_qar100 = 0;dist2_qar200 = 0;dist2_qar400 = 0;dist2_qar800 = 0;

for(i in 1:KK){  
  data1 = qar1(100)
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 = qpsmo(data1,freq,tau, 1)
  
  data2 = qar1(200)
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q2 = qpsmo(data2,freq,tau, 1)
  
  data3 = qar1(400)
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q3 = qpsmo(data3,freq,tau, 1)
  
  data4 = qar1(800)
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q4 = qpsmo(data4,freq,tau, 1)
  
  
  for(j in 1:46){
    q1[,j] = smooth.spline(q1[,j])$y
    q2[,j] = smooth.spline(q2[,j])$y
    q3[,j] = smooth.spline(q3[,j])$y
    q4[,j] = smooth.spline(q4[,j])$y
  }
  q1 = spec.normalize(q1)
  q2 = spec.normalize(q2)
  q3 = spec.normalize(q3)
  q4 = spec.normalize(q4)
  
  q1[q1 < 0] <- 0
  q2[q2 < 0] <- 0
  q3[q3 < 0] <- 0
  q4[q4 < 0] <- 0
  
  dist_qar100 = dist_qar100 + JSD(q1,GT[[5]])/KK
  dist_qar200 = dist_qar200 + JSD(q2,GT[[6]])/KK
  dist_qar400 = dist_qar400 + JSD(q3,GT[[7]])/KK
  dist_qar800 = dist_qar800 + JSD(q4,GT[[8]])/KK
  
  dist2_qar100 = dist2_qar100 + rmse(q1,GT[[5]])/KK
  dist2_qar200 = dist2_qar200 + rmse(q2,GT[[6]])/KK
  dist2_qar400 = dist2_qar400 + rmse(q3,GT[[7]])/KK
  dist2_qar800 = dist2_qar800 + rmse(q4,GT[[8]])/KK
  
  if(i %% 20 == 0){
    print(c('qar',i))}
}


####################################  arfima ##########################################
dist_mix100 = 0;dist_mix200 = 0;dist_mix400 = 0;dist_mix800 = 0;
dist2_mix100 = 0;dist2_mix200 = 0;dist2_mix400 = 0;dist2_mix800 = 0;

for(i in 1:KK){  
  data1 = fracdiff::fracdiff.sim(n = 100, d = 0.2)$series
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 = qpsmo(data1,freq,tau, 1)
  
  data2 = fracdiff::fracdiff.sim(n = 200, d = 0.2)$series
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q2 = qpsmo(data2,freq,tau, 1)
  
  data3 = fracdiff::fracdiff.sim(n = 400, d = 0.2)$series
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q3 = qpsmo(data3,freq,tau, 1)
  
  data4 = fracdiff::fracdiff.sim(n = 800, d = 0.2)$series
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q4 = qpsmo(data4,freq,tau, 1)
  
  
  for(j in 1:46){
    q1[,j] = smooth.spline(q1[,j])$y
    q2[,j] = smooth.spline(q2[,j])$y
    q3[,j] = smooth.spline(q3[,j])$y
    q4[,j] = smooth.spline(q4[,j])$y
  }
  q1 = spec.normalize(q1)
  q2 = spec.normalize(q2)
  q3 = spec.normalize(q3)
  q4 = spec.normalize(q4)
  
  q1[q1 < 0] <- 0
  q2[q2 < 0] <- 0
  q3[q3 < 0] <- 0
  q4[q4 < 0] <- 0
  
  dist_mix100 = dist_mix100 + JSD(q1,GT[[9]])/KK
  dist_mix200 = dist_mix200 + JSD(q2,GT[[10]])/KK
  dist_mix400 = dist_mix400 + JSD(q3,GT[[11]])/KK
  dist_mix800 = dist_mix800 + JSD(q4,GT[[12]])/KK
  
  dist2_mix100 = dist2_mix100 + rmse(q1,GT[[9]])/KK
  dist2_mix200 = dist2_mix200 + rmse(q2,GT[[10]])/KK
  dist2_mix400 = dist2_mix400 + rmse(q3,GT[[11]])/KK
  dist2_mix800 = dist2_mix800 + rmse(q4,GT[[12]])/KK
  
  if(i %% 20 == 0){
    print(c('MIX',i))}
}

###################################  GARCH  ##########################################
dist_garch100 = 0;dist_garch200 = 0;dist_garch400 = 0;dist_garch800 = 0;
dist2_garch100 = 0;dist2_garch200 = 0;dist2_garch400 = 0;dist2_garch800 = 0;

for(i in 1:KK){  
  data1 = as.numeric(garchSim(spec1, n = 100))
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 = qpsmo(data1,freq,tau, 1)
  
  data2 = as.numeric(garchSim(spec1, n = 200))
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q2 = qpsmo(data2,freq,tau, 1)
  
  data3 = as.numeric(garchSim(spec1, n = 400))
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q3 = qpsmo(data3,freq,tau, 1)
  
  data4 = as.numeric(garchSim(spec1, n = 800))
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q4 = qpsmo(data4,freq,tau, 1)
  
  
  for(j in 1:46){
    q1[,j] = smooth.spline(q1[,j])$y
    q2[,j] = smooth.spline(q2[,j])$y
    q3[,j] = smooth.spline(q3[,j])$y
    q4[,j] = smooth.spline(q4[,j])$y
  }
  q1 = spec.normalize(q1)
  q2 = spec.normalize(q2)
  q3 = spec.normalize(q3)
  q4 = spec.normalize(q4)
  
  q1[q1 < 0] <- 0
  q2[q2 < 0] <- 0
  q3[q3 < 0] <- 0
  q4[q4 < 0] <- 0
  
  dist_garch100 = dist_garch100 + JSD(q1,GT[[13]])/KK
  dist_garch200 = dist_garch200 + JSD(q2,GT[[14]])/KK
  dist_garch400 = dist_garch400 + JSD(q3,GT[[15]])/KK
  dist_garch800 = dist_garch800 + JSD(q4,GT[[16]])/KK
  
  dist2_garch100 = dist2_garch100 + rmse(q1,GT[[13]])/KK
  dist2_garch200 = dist2_garch200 + rmse(q2,GT[[14]])/KK
  dist2_garch400 = dist2_garch400 + rmse(q3,GT[[15]])/KK
  dist2_garch800 = dist2_garch800 + rmse(q4,GT[[16]])/KK
  
  if(i %% 20 == 0){
    print(c('garch',i))}
}


resarma = c(dist_arma100,dist_arma200,dist_arma400,dist_arma800,dist2_arma100,dist2_arma200,dist2_arma400,dist2_arma800)
resgarch = c(dist_garch100,dist_garch200,dist_garch400,dist_garch800,dist2_garch100,dist2_garch200,dist2_garch400,dist2_garch800)
resmix = c(dist_mix100,dist_mix200,dist_mix400,dist_mix800,dist2_mix100,dist2_mix200,dist2_mix400,dist2_mix800)
resqar = c(dist_qar100,dist_qar200,dist_qar400,dist_qar800,dist2_qar100,dist2_qar200,dist2_qar400,dist2_qar800)
resmat = cbind(resarma,resqar,resmix,resgarch)

set.panel(1,2)
par(mar=c(3,3,1,1),mgp=c(2,0.5,0)) 
plot(c(100,200,400,800),resarma[1:4],ylim=c(0,0.035),type='l', xlab='', ylab='',
     xaxt='n',yaxt='n',lwd=1.45)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
title(xlab = "n", line = 1.5, ylab = 'RMSE',cex.lab=1) 
lines(c(100,200,400,800),resqar[1:4],col='blue',lwd=1.45)
lines(c(100,200,400,800),resmix[1:4],col='#e6ab02',lwd=1.45)
lines(c(100,200,400,800),resgarch[1:4],col='#1b9e77',lwd=1.45)
points(c(100,200,400,800),resarma[1:4],pch=1)
points(c(100,200,400,800),resqar[1:4],pch=1)
points(c(100,200,400,800),resmix[1:4],pch=1)
points(c(100,200,400,800),resgarch[1:4],pch=1)
legend(x = 'topright', 
       legend = c('ARMA','QAR','ARFIMA','GARCH'),
       lty = 1.45, pch=1,
       col = c('black','blue','#e6ab02','#1b9e77') ,inset = 0.01)

plot(c(100,200,400,800),resarma[5:8],ylim=c(0,0.01),type='l', xlab='', ylab='',
     xaxt='n',yaxt='n',lwd=1.45)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
title(xlab = "n", line = 1.5, ylab = 'JSD',cex.lab=1) 
lines(c(100,200,400,800),resqar[5:8],col='blue',lwd=1.45)
lines(c(100,200,400,800),resmix[5:8],col='#e6ab02',lwd=1.45)
lines(c(100,200,400,800),resgarch[5:8],col='#1b9e77',lwd=1.45)
points(c(100,200,400,800),resarma[5:8],pch=1)
points(c(100,200,400,800),resqar[5:8],pch=1)
points(c(100,200,400,800),resmix[5:8],pch =1)
points(c(100,200,400,800),resgarch[5:8],pch=1)
legend(x = 'topright', 
       legend = c('ARMA','QAR','ARFIMA','GARCH'),
       lty = 1.45, pch=1,
       col = c('black','blue','#e6ab02','#1b9e77') ,inset = 0.01)