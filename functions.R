library(irlba)
library(HDInterval)

############ estimation of pi(all units have the same pi)####################
epsest.func <- function(x,u,sigma){
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  # estimate alternative proportion
  # this implemets the estimator in Jin and Cai
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0, tmax, 0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) { 
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}


#' Estimation of the non-null proportion in Storey 2002
#'
#' Estimates the proportion of non-nulls.
#' 
#' @param pval the corresponding vector of p-values
#' @param lambda the pre-chosen tunning parameter in (0,1)
#' @return a value indicating the estimated non-null proportion
#'
#' @examples
#' x <- rep(c(0,2),c(800,200))+rnorm(1000);
#' pval <- 2*pnorm(-abs(x));
#' epsest.func(pval,0.2);
#'
#' @export
pi.storey <- function(pval,lambda){
  k = as.numeric(sum(pval>lambda))
  pi = 1-k/((1-lambda)*length(pval))
  return(pi)
}

#############functions for multiple testing#########
sc.func<-function(score,score0,q){
  ## USAGE
  # mt.sc(score,score0,q)
  
  ## ARGUMENTS
  # score: original data's score sequence
  # score0: calibration data's score sequence
  # q: the FDR level
  
  ## VALUES
  # nr: the number of rejected hypotheses
  # th: the threshold
  # re: the rejected hypotheses
  # ac: the accepted hypotheses
  # de: the decision rule
  
  m = length(score);
  hps <- rep(0, m);
  
  candidate<-c()
  for (i in 1:m){
    s0 = score0<=score[i] & score>=score0;
    s = score<=score[i] & score<=score0;
    if ((1+sum(s0))/(max(sum(s),1)) <= q){
      candidate<-c(candidate,score[i])
    }
  }
  if (length(candidate)==0){
    threshold<-Inf
    reject<-NULL
    accept<-1:m    
  }else{
    threshold<-max(candidate)
    reject<-which(score<=threshold & score<=score0)
  }
  hps[reject]<-1
  
  y<-list(th=threshold, re=reject, de=hps)
  return (y)
}

bc.func<-function(W,q){
  m=length(W)
  # The decision rule:
  hps <- rep(0, m)
  candidate <- c()
  for (i in 1:m) {
    denom <- (W <= -abs(W[i]))
    numer <- (W >= abs(W[i]))
    if ((1+sum(denom)) / max(sum(numer), 1) <= q) {
      candidate <- c(candidate, abs(W[i]))
    }
  }
  
  if (length(candidate) == 0) {
    threshold <- Inf
    reject <- integer(0)
  } else {
    threshold <- min(candidate)
    reject <- which(W >= threshold)
  }
  
  hps[reject] <- 1
  
  y <- list(th = threshold, re = reject, de = hps)
  return(y)
}

bh.func <- function(pv, q) { 
  # The input
  # pv: the p-values
  # q: the FDR level
  # The output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule
  
  m <- length(pv)
  st.pv <- sort(pv)  # Sorted p-values
  pvi <- st.pv / seq_along(st.pv)  # p-values divided by their rank
  hps <- rep(0, m)  # Initialize decision rule
  
  # Check if all p-values are larger than the threshold
  if (all(pvi > (q/m))) {
    k <- 0
    pk <- 1
    reject <- integer(0)  # No rejection
    accept <- 1:m  # Accept all hypotheses
  } else {
    k <- max(which(pvi <= (q/m)), na.rm = TRUE)
    pk <- st.pv[k]
    reject <- which(pv <= pk)
    accept <- which(pv > pk)
    hps[reject] <- 1  # Mark rejected hypotheses
  }
  
  # Return the result
  y <- list(
    num_rejected = k, 
    threshold = pk, 
    rejected = reject, 
    accepted = accept, 
    decision_rule = hps
  )
  return(y)
}

ad.func<-function(score,score0,q){
  #calculate the conformal p values
  p.ad<-rep(0,m)
  for (i in 1:m){
    p.ad[i] = (1+sum(score[i]>score0))/(m+1)
  }
  y <- bh.func(p.ad,q)
  return(y)
}


ress.func <- function(X, q) {
  # Validate input types
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a matrix or data frame')
  }
  
  # Validate input dimensions
  n <- ncol(X)
  m <- nrow(X)
  stopifnot(n > 1)
  
  # Split the observations into two parts
  n1 <- ceiling(n / 2)
  n2 <- n - n1
  n11 <- sample(1:n, size = n1, replace = F)
  X1 <- X[, n11]
  X2 <- X[, -n11]
  
  if (n >= 4) {
    # Calculate the mean of X1, X2 for each location
    X1.mean <- rowMeans(X1)
    X2.mean <- rowMeans(X2)
    
    # Calculate the standard deviation of X1, X2 for each location
    X1.var <- apply(X1, 1, var)
    X2.var <- apply(X2, 1, var)
  } else if (n == 3) {
    # Calculate the mean of X1, X2 for each location
    X1.mean <- rowMeans(X1)
    X2.mean <- X2
    # Calculate the standard deviation of X1, X2 for each location
    X1.var <- apply(X1, 1, var)
    X2.var <- rep(1, m)
  } else {
    X1.mean <- X1
    X2.mean <- X2
    X1.var <- rep(1, m)
    X2.var <- rep(1, m)
  }
  
  t1 <- sqrt(n1) * X1.mean / sqrt(X1.var)
  t2 <- sqrt(n2) * X2.mean / sqrt(X2.var)
  W <- t1 * t2

  y <- bc.func(W,q)
  return(y)
}

arias.func <- function(X, q) {
  # Validate input types
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a matrix or data frame')
  }
  
  # Validate input dimensions
  n <- ncol(X)
  m <- nrow(X)
  stopifnot(n > 1)
  
  X.mean <- rowMeans(X)
  X.var <- apply(X, 1, var)
  t <- sqrt(n) * X.mean / sqrt(X.var)
  y <- bc.func(t,q)
  return(y)
}
arias1.func <- function(X, q) {
  # Validate input types
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a matrix or data frame')
  }
  
  # Validate input dimensions
  n <- ncol(X)
  m <- nrow(X)
  stopifnot(n > 1)
  
  X.mean <- rowMeans(X)
  X.var <- apply(X, 1, var)
  t <- sqrt(n) * X.mean / sqrt(X.var)
  t0 <- rt(m,n-1)
  t<-sapply(1:m, function(i) sign(t[i]-t0[i])*(max(exp(t[i]),exp(t0[i]))))
  y <- bc.func(t,q)
  return(y)
}
arias2.func <- function(X, q) {
  # Validate input types
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a matrix or data frame')
  }
  
  # Validate input dimensions
  n <- ncol(X)
  m <- nrow(X)
  stopifnot(n > 1)
  
  X.mean <- rowMeans(X)
  X.var <- apply(X, 1, var)
  t <- sqrt(n) * X.mean / sqrt(X.var)
  t0 <- rt(m,n-1)
  
  numerator <- dt(c(t,t0),n-1)
  denominator <- (1-0.1)*dt(c(t,t0),n-1)+0.1*dt(c(t,t0),n-1,sqrt(n));
  Est <- numerator[1:m]/denominator[1:m];
  Est0 <- numerator[(m+1):(2*m)]/denominator[(m+1):(2*m)];
  t <- sapply(1:m, function(i) sign(Est0[i]-Est[i])*max(exp(-Est0[i]),exp(-Est[i])))
  y <- bc.func(t,q)
  return(y)
}
arias3.func <- function(X, q) {
  # Validate input types
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a matrix or data frame')
  }
  
  # Validate input dimensions
  n <- ncol(X)
  m <- nrow(X)
  stopifnot(n > 1)
  
  X.mean <- rowMeans(X)
  X.var <- apply(X, 1, var)
  t <- sqrt(n) * X.mean / sqrt(X.var)
  t0 <- rt(m,n-1)
  
  Est <- 2*pt(-abs(t),n-1);
  Est0 <- 2*pt(-abs(t0),n-1);
  t <- sapply(1:m, function(i) sign(Est0[i]-Est[i])*max(exp(-Est0[i]),exp(-Est[i])))
  
  y <- bc.func(t,q)
  return(y)
}

SENS_oracle <- function(X,alpha,mu,sig){
  #Validate input types.
  if (is.data.frame(X)){
    X.names=names(X)
    X = as.matrix(X,rownames.force=F)
  } else if (is.matrix(X)){
    X.names=colnames(X)
  } else{
    stop('Input X must be a matrix or data frame')}
  
  #Validate input dimensions.
  n <- ncol(X); m <- nrow(X);
  stopifnot(n>1)
  
  #Split the observations into two parts
  n1 <- ceiling(n/2); n2=n-n1;
  n11 <- sample(1:n,size=n1,replace=F); 
  X1 <- X[,n11]; X2 <- X[,-n11];
  
  if(n>=4){
    #Calculate the mean of X1, X2 for each location
    X1.mean <- rowMeans(X1);
    X2.mean <- rowMeans(X2);
    
    V <- sqrt(n1*n2/n)*(X1.mean+X2.mean);
    V0 <- sqrt(n1*n2/n)*(X1.mean-X2.mean);
    #Calculate the standard deviation of X1, X2 for each location
    X1.var <- apply(X1,1,var);
    X2.var <- apply(X2,1,var);
    #Calculate the main statistics used later
    S <- sqrt(((n1-1)*X1.var+(n2-1)*X2.var)/(n-2));
  }else if(n==3){
    #Calculate the mean of X1, X2 for each location
    X1.mean <- rowMeans(X1);
    X2.mean <- X2;
    
    V <- sqrt(n1*n2/n)*(X1.mean+X2.mean);
    V0 <- sqrt(n1*n2/n)*(X1.mean-X2.mean);
    #Calculate the standard deviation of X1, X2 for each location
    X1.var <- apply(X1,1,var);
    #Calculate the main statistics used later
    S <- X1.var;
  }else{
    V <- sqrt(n1*n2/n)*(X1+X2);
    V0 <- sqrt(n1*n2/n)*(X1-X2);
    S <-rep(1,m);
  }
  t <- V/S; t0 <- V0/S;
  
  sens.numerator <- dt(c(t,t0),n-2)
  sens.denominator <- (1-0.1)*dt(c(t,t0),n-2)+0.1*dt(c(t,t0),n-2,2*sqrt(n1*n2/n)*mu/sig);
  sens.Est <- sens.numerator[1:m]/sens.denominator[1:m];
  sens.Est0 <- sens.numerator[(m+1):(2*m)]/sens.denominator[(m+1):(2*m)];
  t<-sapply(1:m, function(i) sign(sens.Est0[i]-sens.Est[i])*(max(exp(-sens.Est[i]),exp(-sens.Est0[i]))))
  y = bc.func(t,alpha)
  return(y)
}


############functions for estimation of empirical null#############
EstNull.func<-function (x,gamma=0.1){
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation
  
  n = length(x)
  t = c(1:1000)/200
  
  gan    = n^(-gamma)
  that   = 0 
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)
  
  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]
  
  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat) 
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)
  
  return(musigma=list(mu=uhat,s=shat))
}









