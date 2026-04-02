source('fn.r')
library(quantreg)
dat <- readMat('IR007_0.mat')
dat = dat$X105.DE.time[,1]

d1= dat[1:2000]
n = length(d1)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.5]
sqp1 = spec.normalize( qpsmo(d1, freq, 0.8, sd(d)))
pg1 = spec.normalize( spec.pgram(d1,plot=F)$spec)

d2= dat[2001:4000]
n = length(d2)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.5]
sqp2 = spec.normalize( qpsmo(d2, freq, 0.8, sd(d)))
pg2 = spec.normalize( spec.pgram(d2,plot=F)$spec)

d3= dat[4001:6000]
n = length(d3)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.5]
sqp3 = spec.normalize( qpsmo(d3, freq, 0.8, sd(d)))
pg3 = spec.normalize( spec.pgram(d3,plot=F)$spec)

d4= dat[6001:8000]
n = length(d4)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.5]
sqp4 = spec.normalize( qpsmo(d4, freq, 0.8, sd(d)))
pg4 = spec.normalize( spec.pgram(d4,plot=F)$spec)


##########plot
par(mar=c(3,3,2,1),mgp=c(2,0.5,0)) 
layout_matrix <- matrix(c(1, 1, 2, 3, 4, 4,5,6,7,7,8,9,10,10,11,12), nrow =4, byrow = TRUE)
layout(mat = layout_matrix, widths = c(1, 1 ,1,1))

plot(d1, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("Segment 1"),line=1)
title(xlab = "Time index", line = 1.5,cex.lab=1.2) 
grid()
lines(-0.35* cos( 2*pi*freq[27]*((1:2000)) )+0.5 ,col='blue')
lines( 0.35* cos( 2*pi*freq[27]*((1:2000)) )-0.5 ,col='blue')

plot(freq, sqp1, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("SQP"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()

plot(freq, pg1[2:1000], type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("PG"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()


plot(d2, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("Segment 2"),line=1)
title(xlab = "Time index", line = 1.5,cex.lab=1.2) 
grid()
lines(-0.35* cos( 2*pi*freq[27]*((1:2000)) )+0.5 ,col='blue')
lines( 0.35* cos( 2*pi*freq[27]*((1:2000)) )-0.5 ,col='blue')

plot(freq, sqp2, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("SQP"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()

plot(freq, pg2[2:1000], type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("PG"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()



plot(d3, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("Segment 3"),line=1)
title(xlab = "Time index", line = 1.5,cex.lab=1.2) 
grid()
lines(-0.35* cos( 2*pi*freq[27]*((1:2000)) )+0.5 ,col='blue')
lines( 0.35* cos( 2*pi*freq[27]*((1:2000)) )-0.5 ,col='blue')

plot(freq, sqp3, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("SQP"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()

plot(freq, pg3[2:1000], type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("PG"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()


plot(d4, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("Segment 4"),line=1)
title(xlab = "Time index", line = 1.5,cex.lab=1.2) 
grid()
lines(-0.35* cos( 2*pi*freq[27]*((1:2000)) )+0.5 ,col='blue')
lines( 0.35* cos( 2*pi*freq[27]*((1:2000)) )-0.5 ,col='blue')

plot(freq, sqp4, type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("SQP"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()

plot(freq, pg4[2:1000], type='l', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(main = expression("PG"),line=1)
title(xlab = "Freq", line = 1.5,cex.lab=1.2) 
grid()


#######################
dat <- readMat('IR007_0.mat')
dat = dat$X105.DE.time[,1]

tau = c(0.1, 0.5, 0.9)
win_len <- 2000
hop <- 1000
n <- length(dat)
starts <- seq(1, n - win_len + 1, by = hop)
win_index <- lapply(starts, function(s) c(start = s, end = s + win_len - 1))

qp1 = list()
pg1 = list()
for(i in 1:length(win_index)){
  data1 = dat[win_index[[i]][1]:win_index[[i]][2]]
  n = length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 <- spec.normalize( qpsmo(data1, freq, tau, sd(data1)))
  qp1[[i]] = q1
  pg1[[i]] = spec.pgram(data1,plot=F)$spec[1:999]
  print(i)
}


corpg1 = rep(length(pg1),0)
corqp11 = rep(length(pg1),0)
corqp21 = rep(length(pg1),0)
corqp31 = rep(length(pg1),0)


for(i in 1:(length(pg)-1)){
  corpg1[i] = cor(pg1[[i]],pg1[[i+1]])
  corqp11[i] = cor(qp1[[i]][,1],qp1[[i+1]][,1])
  corqp21[i] = cor(qp1[[i]][,2],qp1[[i+1]][,2])
  corqp31[i] = cor(qp1[[i]][,3],qp1[[i+1]][,3])
}

set.panel(1,2)
par(mar=c(3,3,1,1),mgp=c(2,0.5,0)) 
plot(corpg1,type='l', xlab='', ylab='',xaxt='n',yaxt='n',lwd=1.35,ylim=c(0.75,1))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
title(xlab = "Window index", line = 1.5, ylab = 'Correlation',cex.lab=1.2) 
lines(corqp11,col='blue',lwd=1.35)
lines(corqp21,col='#e6ab02',lwd=1.35)
lines(corqp31,col='#1b9e77',lwd=1.35)
legend(x = 'bottomright', 
       legend = c(expression('PG'),
                  expression( 'SQP('~tau==0.1~')' ),expression( 'SQP('~tau==0.5~')' ),
                  expression( 'SQP('~tau==0.9~')' )),
       lty = c(1,1,1,1),
       col = c( 'black','blue','#e6ab02','#1b9e77'),inset = 0.01,ncol=2)

boxplot(
  cbind(corpg1, corqp11, corqp21, corqp31),
  xaxt = "n"
)
axis(1, at = 1:4, labels = c('PG', 'SQP (0.1)', 'SQP (0.5)', 'SQP (0.9)'), tck = -0.02)
axis(2, tck = -0.02)



