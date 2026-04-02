rho_tc <- function(u, tau, c) {
  ((2 * tau - 1) * u + sqrt(c^2 + u^2)) / 2
}

rho_t <- function(u, tau) {
  ifelse(u >= 0, tau * u, (tau - 1) * u)
}

smoqr <- function(y, x1, x2, c = NULL, tau, c_kappa = 1) {
  if (length(y) != length(x1) || length(y) != length(x2)) {
    stop("y, x1, x2 must have the same length.")
  }
  if (!is.null(c) && (!is.numeric(c) || length(c) != 1 || c <= 0)) {
    stop("c must be NULL or a positive number.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0 || tau >= 1) {
    stop("tau must be a number in (0, 1).")
  }
  if (!is.numeric(c_kappa) || length(c_kappa) != 1 || c_kappa <= 0) {
    stop("c_kappa must be a positive number.")
  }
  
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  dat <- dat[complete.cases(dat), ]
  
  y <- dat$y
  X <- cbind(1, dat$x1, dat$x2)
  n <- length(y)
  
  # initial value: OLS with intercept
  beta_init <- tryCatch(
    as.vector(lm.fit(x = X, y = y)$coefficients),
    error = function(e) c(0, 0, 0)
  )
  beta_init[!is.finite(beta_init)] <- 0
  
  # auto select c if c is NULL
  if (is.null(c)) {
    r0 <- as.vector(y - X %*% beta_init)
    
    # robust scale from residuals
    s0 <- mad(r0, center = median(r0), constant = 1.4826)
    
    # fallbacks
    if (!is.finite(s0) || s0 <= 0) s0 <- sd(r0)
    if (!is.finite(s0) || s0 <= 0) s0 <- mad(y, center = median(y), constant = 1.4826)
    if (!is.finite(s0) || s0 <= 0) s0 <- sd(y)
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
    
    c_use <- c_kappa * s0 * n^(-1/3)
  } else {
    c_use <- c
  }
  
  obj_fun <- function(beta) {
    r <- as.vector(y - X %*% beta)
    sum(((2 * tau - 1) * r + sqrt(c_use^2 + r^2)) / 2)
  }
  
  grad_fun <- function(beta) {
    r <- as.vector(y - X %*% beta)
    psi <- ((2 * tau - 1) + r / sqrt(c_use^2 + r^2)) / 2
    -as.vector(crossprod(X, psi))
  }
  
  fit <- optim(
    par = beta_init,
    fn = obj_fun,
    gr = grad_fun,
    method = "BFGS"
  )
  
  beta_hat <- fit$par
  names(beta_hat) <- c("intercept", "beta_x1", "beta_x2")
  
  return(list(
    coefficients = beta_hat[2:3],
    intercept = beta_hat[1],
    c = c_use,
    objective = fit$value,
    convergence = fit$convergence
  ))
}


qpsmo <- function(y, f, tau, cc) {
  ns <- length(y)
  nf <- length(f)
  ntau <- length(tau)
  tt <- 1:ns
  
  out <- matrix(NA_real_, nrow = nf, ncol = ntau)
  
  for (k in 1:nf) {
    ff <- f[k]
    x1 <- cos(2 * pi * ff * tt)
    x2 <- sin(2 * pi * ff * tt)
    
    coef_mat <- matrix(NA_real_, nrow = 2, ncol = ntau)
    for (i in 1:ntau) {
      res = smoqr(y, x1, x2, cc, tau[i])
      coef_mat[, i] <- res$coefficients
    }
    
    out[k, ] <- colSums(coef_mat^2) * ns / 4
  }
  
  out[out < 0] <- 0
  if (ntau == 1) out <- c(out)
  out
}

qpsmoc <- function(y, f, tau) {
  ns <- length(y)
  nf <- length(f)
  ntau <- length(tau)
  tt <- 1:ns
  
  out <- matrix(NA_real_, nrow = nf, ncol = ntau)
  
  ff <- f[1]
  x1 <- cos(2 * pi * ff * tt)
  x2 <- sin(2 * pi * ff * tt)
  cc = smoqr(y, x1, x2, NULL, tau[1])$c
  print
  
  for (k in 1:nf) {
    ff <- f[k]
    x1 <- cos(2 * pi * ff * tt)
    x2 <- sin(2 * pi * ff * tt)
    
    coef_mat <- matrix(NA_real_, nrow = 2, ncol = ntau)
    for (i in 1:ntau) {
      res = smoqr(y, x1, x2, cc, tau[i])
      coef_mat[, i] <- res$coefficients
    }
    
    out[k, ] <- colSums(coef_mat^2) * ns / 4
  }
  
  out[out < 0] <- 0
  if (ntau == 1) out <- c(out)
  return(list(out=out,c=cc))
}


expectile_peri_ls <- function(y, f, tau, n.cores=1) {
  ep.parallel1<-function(yy,ff,tau,tt) {
    ntau = length(tau)
    coef = matrix(0, 2, ntau)
    for(i in 1:ntau){
      taui = tau[i]
      x1 = cos(2*pi*ff*tt)
      x2 = sin(2*pi*ff*tt)
      o = rep(1,length(x1))
      X = cbind(o,x1,x2)
      fit <- expectreg.ls(yy ~ X[,2]+ X[,3], mstop=200, expectiles=taui, quietly=TRUE)
      coef[,i] = c(fit$coefficients[[1]],fit$coefficients[[2]])
    }
    return(coef)
  }
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  yy=y
  tmp=list()
  for(i in 1:length(f)){
    tmp[[i]] = ep.parallel1(yy,f[i],tau,tt)
  }
  
  out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
  out<-matrix(unlist(out),ncol=ntau,byrow=T)
  out<-out*ns/4
  
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
  
}


lap.spec.new2 <- function(y, f, tau=0.5, intercept=F, type=1, weights=NULL, method="fn", n.cores=1) {
  # type 1: squared L2 norm of coefficients
  # type 2: cost difference
  
  rsoid2 <- function(n, f, a, b) {
    p <- length(f)
    tt <- c(1:n)
    one <- rep(1, n)
    tmp <- (one %o% a) * cos(2*pi*(tt %o% f)) + (one %o% b) * sin(2*pi*(tt %o% f))
    tmp <- apply(tmp,1,sum)
    tmp
  }
  
  
  lap.cost<-function(y,tau=0.5,weights=NULL) {
    # cost function of quantile regression
    ns<-length(y)
    if(is.null(weights)) weights<-rep(1,ns)
    tmp<-tau*y
    sel<-which(y < 0)
    if(length(sel) > 0) tmp[sel]<-(tau-1)*y[sel]
    sum(tmp*weights,na.rm=T)
  }
  
  
  qh.parallel<-function(yy,ff,tau,tt,ns,weights,type) {
    # parallel computation of quantile harmonic regression
    # for a single value of tau (and single value of ff)
    # used when intercept = F
    if(ff==0.5) {
      fit<-try(rq(yy ~ 0+cos(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-c(0,0)
      } else {
        fit$coefficients<-c(fit$coefficients,0)
      }
    }
    if(ff==0) {
      fit<-NULL
      fit$coefficients<-c(0,0)
    }
    if(ff != 0.5 & ff>0) {
      fit<-try(rq(yy ~ 0+cos(2*pi*ff*tt)+sin(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-c(0,0)
      }
    }
    fit$residuals<-yy-rsoid2(ns,ff,fit$coefficients[1],fit$coefficients[2])
    if(type==1) tmp.coef<-fit$coefficients
    if(type==2) tmp.cost<-lap.cost(fit$residuals,tau=tau,weights=weights)
    rm(fit)
    if(type==1) return(tmp.coef)
    if(type==2) return(tmp.cost)
  }
  
  
  qh.parallel2<-function(yy,ff,tau,tt,ns,weights,type) {
    # parallel computation of quantile harmonic regression
    # for a vector of tau (and single value of ff)
    # used when intercept = T
    if(ff==0.5) {
      fit<-try(rq(yy ~ cos(2*pi*ff*tt),method=method,tau=tau,weights=weights),silient=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
      } else {
        fit$coefficients<-rbind(fit$coefficients,rep(0,length(tau)))
      }
    }
    if(ff==0) {
      fit<-NULL
      fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
    }
    if(ff != 0.5 & ff > 0) {
      fit<-try(rq(yy ~ cos(2*pi*ff*tt)+sin(2*pi*ff*tt),method=method,tau=tau,weights=weights),silent=T)
      if(length(fit)==1) {
        fit<-NULL
        fit$coefficients<-rbind(quantile(yy,probs=tau),rep(0,length(tau)),rep(0,length(tau)))
      }
    }
    fit$coefficients<-matrix(fit$coefficients,ncol=length(tau))
    if(type==1) tmp.coef<-matrix(fit$coefficients[-1,],ncol=length(tau))
    if(type==2) {
      if(length(tau)==1) fit$coefficients<-matrix(fit$coefficients,ncol=1)
      tmp.cost<-rep(NA,length(tau))
      for(i.tau in c(1:length(tau))) {
        tmp.resid<-yy-fit$coefficients[1,i.tau]-rsoid2(ns,ff,fit$coefficients[2,i.tau],fit$coefficients[3,i.tau])
        tmp.cost[i.tau]<-lap.cost(tmp.resid,tau=tau[i.tau],weights=weights)
      }
      rm(tmp.resid)
    }
    rm(fit)
    if(type==1) return(tmp.coef)
    if(type==2) return(tmp.cost)
  }
  
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  if(n.cores>1) {
    library(foreach) 
    library(doParallel) # libraries required for parallel computing
    classify.cl <- makeCluster(n.cores)
    registerDoParallel(classify.cl) 
  }
  
  
  if(is.null(weights)) weights<-rep(1,ns)
  
  
  yy<-y
  if(intercept) {
    
    if(type==1) coef<-array(NA,dim=c(ntau,2,nf))
    if(type==2) { 
      cost<-matrix(NA,ntau,nf)
      fit<-rq(yy ~ 1,method=method,tau=tau,weights=weights)
      if(ntau==1) fit$coefficients<-matrix(fit$coefficients,ncol=1)
      cost0<-rep(NA,ntau)
      for(i.tau in c(1:ntau)) cost0[i.tau]<-lap.cost(yy-fit$coefficients[,i.tau],tau=tau[i.tau],weights=weights)
    }
    
    if(n.cores>1) {
      tmp<-foreach(i=1:nf,.packages="quantreg") %dopar% { qh.parallel2(yy,f[i],tau=tau,tt,ns,weights,type) }
    } else {
      library(foreach)
      tmp<-foreach(i=1:nf,.packages="quantreg") %do% { qh.parallel2(yy,f[i],tau=tau,tt,ns,weights,type) }
    }
    if(type==1) {
      out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
      out<-matrix(unlist(out),ncol=ntau,byrow=T)
      out<-out*ns/4
    }
    if(type==2) {
      out<-matrix(unlist(tmp),ncol=ntau,byrow=T)
      out<-matrix(rep(cost0,nf),ncol=ntau,byrow=T)-out
    }
    
  } else {
    
    out<-NULL
    for(i.tau in c(1:ntau)) {
      if(type==1) coef<-matrix(NA,ncol=2,nrow=nf)
      if(type==2) yy<-y-quantile(y,probs=tau[i.tau],na.rm=T)
      cost0<-lap.cost(yy,tau=tau[i.tau],weights=weights)
      if(n.cores>1) {
        tmp<-foreach(i=1:nf,.packages="quantreg") %dopar% { qh.parallel(yy,f[i],tau=tau[i.tau],tt,ns,weights,type) }
      } else {
        library(foreach)
        tmp<-foreach(i=1:nf,.packages="quantreg") %do% { qh.parallel(yy,f[i],tau=tau[i.tau],tt,ns,weights,type) }
      }
      if(type==1) {
        coef<-matrix(unlist(tmp),ncol=2,byrow=T)
        out<-cbind(out,apply(coef^2,1,sum)*ns/4)
      }
      if(type==2) {
        cost<-c(unlist(tmp))
        out<-cbind(out,cost0-cost)
      }
    }
    
  }
  
  if(n.cores>1) stopCluster(classify.cl)
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
  
}


pfish = function(x){
  q = length(x)
  g = q * max(x)/sum(x)
  js = 0:q
  s = 1-js*g/q
  for(i in 1:(q+1)){
    if (s[i] <0) s[i] = 0
  }
  pp = 1 - sum((-1)^js * choose(q,js) * s^(q-1))
  return(pp)
}


spec.normalize<-function(qper) {
  if(is.matrix(qper)) {
    return( apply(qper,2,FUN=function(x) { x/sum(x) }) )
  } else {
    return( qper/sum(qper) )
  }
}


mse_mat <- function(A, B, na.rm = FALSE) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (!all(dim(A) == dim(B))) {
    stop("A and B must have the same dimensions.")
  }
  diff <- A - B
  sqrt( mean(diff^2, na.rm = TRUE))
}



type1 = function(yt,out){
  yt[floor(runif(1,length(yt),n=1))] = out
  return(yt)}

type2 = function(yt,out,len){
  ind = floor(runif(1,length(yt)-len,n=1))
  yt[ind:(ind+len-1)] = yt[ind:(ind+len-1)]+ out
  return(yt)}


type3 = function(yt,out,len){
  ind = floor(runif(1,length(yt)-len,n=1))
  yt[ind:(ind+len-1)] = yt[ind:(ind+len-1)]*out
  return(yt)}


# type3 = function(xt, amp){
#   x = xt[1:(length(xt)/2)]
#   len = length(x)
# 
#   c = floor(0.031*len):(floor(0.27*len))/(0.15*len)
#   ga = gamma(c)-gamma(1.8)
#   ga_1 = 1.5*amp*ga[length(ga):1]
#   ga_2 = -amp*ga
#   ga_3 = c(ga_1,ga_2)
#   inter = seq(from=ga_3[floor(length(ga_3)/2)], to = ga_3[floor(length(ga_3)/2)+1],
#               by =  (ga_3[floor(length(ga_3)/2)+1] - ga_3[floor(length(ga_3)/2)]) / floor(0.032*len))
# 
#   gam = c()
#   gam[1:length(ga_1)] = ga_1
#   gam[(length(ga_1)+1):(length(ga_1)+length(inter)-2)] = inter[2:(length(inter)-1)]
#   gam[(length(ga_1)+length(inter)-1):((length(ga_1)+length(inter))+length(ga_2)-2)] = ga_2
# 
# 
#   x[(floor(0.25*len)+1):(floor(0.25*len)+length(gam))] =
#     x[(floor(0.25*len)+1):(floor(0.25*len)+length(gam))] + gam
#   xt[1:(length(xt)/2)] = x
#   return(xt)
# }



qar1 <- function(n) {
  u <- runif(n)
  y <- numeric(n)
  for (t in 2:n) {
    y[t] <- 0.1 * qnorm(u[t]) + 1.9 * (u[t] - 0.5) * y[t - 1]
  }
  y
}


rmse <- function(y1, y2) {
  sqrt(mean((y1 - y2)^2))
}



JSD <- function(sqp1, sqp2, base = 2, eps = 1e-12) {
  # check input
  if (!is.matrix(sqp1) || !is.matrix(sqp2)) {
    stop("sqp1 and sqp2 must both be matrices.")
  }
  if (!all(dim(sqp1) == dim(sqp2))) {
    stop("sqp1 and sqp2 must have the same dimensions.")
  }
  if (any(!is.finite(sqp1)) || any(!is.finite(sqp2))) {
    stop("sqp1 and sqp2 must contain only finite values.")
  }
  if (any(sqp1 < 0) || any(sqp2 < 0)) {
    stop("Jensen-Shannon divergence requires nonnegative entries.")
  }
  
  nalpha <- ncol(sqp1)
  
  # KL divergence for probability vectors
  KL <- function(p, q, base = 2) {
    idx <- p > 0
    sum(p[idx] * (log(p[idx] / q[idx]) / log(base)))
  }
  
  jsd_each <- numeric(nalpha)
  
  for (j in 1:nalpha) {
    p <- sqp1[, j]
    q <- sqp2[, j]
    
    # avoid zero-sum columns
    if (sum(p) <= 0 || sum(q) <= 0) {
      stop(sprintf("Column %d has nonpositive total mass.", j))
    }
    
    # normalize each column into a probability vector
    p <- p / sum(p)
    q <- q / sum(q)
    
    # numerical stabilization
    p <- p + eps
    q <- q + eps
    p <- p / sum(p)
    q <- q / sum(q)
    
    m <- 0.5 * (p + q)
    
    jsd_each[j] <- 0.5 * KL(p, m, base = base) + 0.5 * KL(q, m, base = base)
  }
  
return(mean(jsd_each))
}