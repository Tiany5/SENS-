#' The SENS procedure for controlling the false discovery rate
#'
#' @docType package
#' @name CEN
library(irlba);
library(Matrix);
library(iterators);
library(REBayes);
library(ggpubr);
library(rmutil);
library(cubature);
library(Rmosek);
library(densratio);
library(mkde);

#' SENS
#'
#'This function runs the SENS procedure, constructing the calibrated empirical null and the test data, and constructing auxiliary statistics under heteroskedasticity, choosing the cutoff and selectiing the locations.
#'
#' @param X the matrix or data frame of observation
#' @param alpha targeted FDR (false discovery rate) level
#' @param option Gasssian case for the null distribution or General case for the null distribution
#' #'
#' @return A list containing the following components: 
#' \item{de}{decision for each location (0 or 1)}
#' \item{th}{threshold for SENS procedure}
#' 
#' @examples
#' X <- matrix(rep(c(0,3),c(800,200))+rnorm(1000),ncol=5,nrow=200);
#' SENS(X,0.05,'Gaussian');
#'
#' @importFrom density
#' @importFrom stats density dt pnorm pt qnorm var rnorm
#' @export 
SENS2 <- function(X,Y,alpha,option=c('Gaussian','General')){
  #Validate input types.
  if (is.data.frame(X)){
    X.names=names(X)
    X = as.matrix(X,rownames.force=F)
  } else if (is.matrix(X)){
    X.names=colnames(X)
  } else{
    stop('Input X must be a matrix or data frame')}
  
  #Validate input dimensions.
  nx <- ncol(X); mx <- nrow(X);
  stopifnot(nx>1)
  
  #Split the observations into two parts
  nx1 <- ceiling(nx/2); nx2=nx-nx1;
  nx11 <- sample(1:nx,size=nx1,replace=F); 
  X1 <- X[,nx11]; X2 <- X[,-nx11];
  
  #Calculate the mean of X1, X2 for each location
  X1.mean <- rowMeans(X1);
  X2.mean <- rowMeans(X2);
    
  Vx <- (X1.mean+X2.mean);
  Vx0 <-(X1.mean -X2.mean);
  #Calculate the standard deviation of X1, X2 for each location
  X1.var <- apply(X1,1,var);
  X2.var <- apply(X2,1,var);
  #Calculate the main statistics used later
  Sx <- sqrt(nx/nx1*nx2)*sqrt(((nx1-1)*X1.var+(nx2-1)*X2.var)/(nx-2));
  
  #Validate input types.
  if (is.data.frame(X)){
    X.names=names(X)
    X = as.matrix(X,rownames.force=F)
  } else if (is.matrix(X)){
    X.names=colnames(X)
  } else{
    stop('Input X must be a matrix or data frame')}
  
  #Validate input dimensions.
  ny <- ncol(Y); my <- nrow(Y);
  stopifnot(ny>1)
  
  #Split the observations into two parts
  ny1 <- ceiling(ny/2); ny2=ny-ny1;
  ny11 <- sample(1:ny,size=ny1,replace=F); 
  Y1 <- Y[,ny11]; Y2 <- Y[,-ny11];
  
  #Calculate the mean of Y1, Y2 for each location
  Y1.mean <- rowMeans(Y1);
  Y2.mean <- rowMeans(Y2);
  
  Vy <- (Y1.mean+Y2.mean);
  Vy0 <-(Y1.mean -Y2.mean);
  #Calculate the standard deviation of Y1, Y2 for each location
  Y1.var <- apply(Y1,1,var);
  Y2.var <- apply(Y2,1,var);
  #Calculate the main statistics used later
  Sy <- sqrt(ny/ny1*ny2)*sqrt(((ny1-1)*Y1.var+(ny2-1)*Y2.var)/(ny-2));
  
  V<-Vx-Vy;V0<-Vx0-Vy0
  S<-sqrt(Sx^2+Sy^2)
  t <- V/S; t0 <- V0/S;
  
  # Use t and t0 as the test and calibration
  # Calculate the estimated SENS statistics based on density estimation for denominator
  t <- inverseCDF(pt(t,df=nx+ny-4),pnorm); t0 <- inverseCDF(pt(t0,df=nx+ny-4),pnorm);
  
  if (option=='Gaussian'){
    Em<-EstNull.func(c(t,t0))
    density_values0<-dnorm(c(t,t0),Em$mu,Em$s)
  }else{
    t00 <- ifelse(abs(t) < abs(t0), t, t0)
    density_values0 <- numeric(2*m)
    bw0<-density(c(t00,-t00))$bw
    for (i in 1:(2*m)){
      density_values0[i]<-sum(dnorm((c(t,t0)[i]-c(t00,-t00))/bw0))/(2*m*bw0)
    }
  }
  
  sens.numerator <- density_values0
  sens.denominator <- approxfun(density(c(t,t0))$x,density(c(t,t0))$y)(c(t,t0));
  sens.Est <- sens.numerator[1:m]/sens.denominator[1:m];
  sens.Est0 <- sens.numerator[(m+1):(2*m)]/sens.denominator[(m+1):(2*m)];
  t<-sapply(1:m, function(i) sign(sens.Est0[i]-sens.Est[i])*(max(exp(-sens.Est[i]),exp(-sens.Est0[i]))))
  y = bc.func(t,alpha)
  return(y)
}


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



#' Conformal inference for FDR control
#'
#' Decision rule based on conformal inference.
#' 
#' @param score original data's score sequence
#' @param score0 calibration data's score sequence
#' @param alpha targeted FDR (false discovery rate) level
#' @return A list containing the following components: 
#' \item{de}{decision for each location (0 or 1)}
#' \item{th}{threshold for the decision rule}
#'
#' @examples
#' score <- 2*pnorm(-abs(rep(c(0,3),c(800,200))+rnorm(1000)))
#' score0 <- 2*pnorm(-abs(rnorm(1000)))
#' alpha <- 0.05
#' sc.func(score,score0,alpha);
#'
#' @export
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


EstNull.func<-function (x,gamma=0.1){
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation
  
  n = length(x)
  t = c(1:10000)/200
  
  gan    = n^(-gamma)
  that   = 0 
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,10000)
  phiminus  = rep(1,10000)
  dphiplus  = rep(1,10000)
  dphiminus = rep(1,10000)
  phi       = rep(1,10000)
  dphi      = rep(1,10000)
  
  for (i in 1:10000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:10000)[(phi - gan) <= 0])
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



