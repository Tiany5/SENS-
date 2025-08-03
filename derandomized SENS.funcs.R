#' The CEN procedure for controlling the false discovery rate
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

#' CENS
#'
#'This function runs the CEN procedure, constructing the calibrated empirical null and the test data, and constructing auxiliary statistics under heteroskedasticity, choosing the cutoff and selectiing the locations.
#'
#' @param X the matrix or data frame of observation
#' @param alpha targeted FDR (false discovery rate) level
#' @param length_out the number of the grid points of mu
#' @param bases the dimension of the space we construct for the estimate of mu
#' @param option homoskedastic case for FDR control or heteroskedastic for FDR control
#' #'
#' @return A list containing the following components: 
#' \item{de}{decision for each location (0 or 1)}
#' \item{th}{threshold for CARS procedure}
#' 
#' @examples
#' X <- matrix(rep(c(0,3),c(800,200))+rnorm(1000),ncol=5,nrow=200);
#' CEN(X,0.05,50,10,'homoskedastic','working model');
#'
#' @importFrom density
#' @importFrom stats density dt pnorm pt qnorm var rnorm
#' @export 
derandomized_SENS <- function(X,alpha,N,pro,option=c('Gaussian','General')){
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
  e <- matrix(0,nrow = m,ncol = N)
  for(k in 1:N){
  set.seed(k)
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
  t[which(is.na(t))]<-0; t0[which(is.na(t0))]<-0
  
  # Use t and t0 as the test and calibration
  # Calculate the estimated SENS statistics based on density estimation for denominator
  if(n>=3){
    t <- inverseCDF(pt(t,df=n-2),pnorm); t0 <- inverseCDF(pt(t0,df=n-2),pnorm);}
  
  if (option=='Gaussian'){
    Em<-EstNull.func(c(t,t0))
    density_values0<-dnorm(c(t,t0),Em$mu,Em$s)
  }else{
    t00<-c()
    for(i in 1:m){
      if(abs(t[i])>=abs(t0[i])){
        t00<-c(t00,t0[i])
      }else{
        t00<-c(t00,t[i])
      }
    }
    density_values0 <- numeric(2*m)
    bw0<-density(t00)$bw
    for (i in 1:(2*m)){
      density_values0[i]<-sum(dnorm((c(t,t0)[i]-t00)/bw0))/(m*bw0)
    }
  }
  
  sens.numerator <- density_values0
  sens.denominator <- approxfun(density(c(t,t0))$x,density(c(t,t0))$y)(c(t,t0));
  sens.Est <- sens.numerator[1:m]/sens.denominator[1:m];
  sens.Est0 <- sens.numerator[(m+1):(2*m)]/sens.denominator[(m+1):(2*m)];
  y = sc.func(sens.Est,sens.Est0,alpha*pro)
  e[,k]=m*y$de/(1+sum((sens.Est0<=y$th)*(sens.Est0<=sens.Est)))
  }
  e = rowMeans(e)
  y = eBH.func(e,alpha)
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
  for (i in 1:m){
    s0 = score0<=score0[i] & score>=score0;
    s = score<=score0[i] & score<=score0;
    if ((1+sum(s0))/(max(sum(s),1)) <= q){
      candidate<-c(candidate,score0[i])
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



#' eBH for FDR control
#'
#' Decision rule based on eBH.
#' 
#' @param ev e-values
#' @param alpha targeted FDR (false discovery rate) level
#' @return A list containing the following components: 
#' \item{de}{decision for each location (0 or 1)}
#' \item{th}{threshold for the decision rule}
#' \item{nr}{the number of hypothesis to be rejected}
#' \item{re}{the index of rejected hypotheses}
#' \item{ac}{the index of accepted hypotheses}
#'
#' @examples
#' ev <- runif(1000,0,2)
#' alpha <- 0.05
#' eBh.func(ev,alpha);
#'
#' @export
eBH.func <- function(ev,alpha){
  st.ev<-sort(ev,decreasing = T)
  evi<-st.ev*(1:m)/m
  hps<-rep(0, m)
  if (max(evi>=(1/alpha))==0){
    k<-0
    evk<-1
    reject<-NULL
    accept<-1:m
  }else{
    k<-max(which(evi>=(1/alpha)))
    evk<-st.ev[k]
    reject<-which(ev>=evk)
    accept<-which(ev<evk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=evk, re=reject, ac=accept, de=hps)
  return (y)
}
