#' @title Main function for fast estimation CUB models 
#' @description Main function to estimate and validate a CUB model for explaining uncertainty 
#' and feeling for given ratings, with or without covariates, on the basis of Louis' identity for the information matrix and the derived accelerated estimation.
#' @aliases fastCUB
#' @usage fastCUB(Formula, data, ...)
#' @param Formula Object of class Formula, with two matrix for the right-hand side for uncertainty and feeling covariates.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param ... Additional arguments to be passed for the specification of the model, including covariates matrices Y, W, X for 
#' #' for uncertainty, feeling and shelter, respectively.
#' @return An object of the class "fastCUB": returns a list containing the following results: 
#' \item{estimates}{Maximum likelihood estimates of parameters}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' @import Formula
#' @keywords package
#' @export fastCUB
#' @details This is the main function for CUB models, which calls for the corresponding functions whenever 
#' covariates are specified. It performs maximum likelihood estimation via the E-M algorithm 
#' for CUB models and extensions based on the Louis'identity for the observed information matrix.
#' @seealso \code{\link{probcub00}}, \code{\link{probcubp0}}, \code{\link{probcub0q}}, \code{\link{probcubpq}},


fastCUB<-function(Formula,data,...){
  
  call <- match.call()
  if (missing(data)) {
    data <- environment(Formula)
  }
  
  mf <- model.frame(Formula, data, na.action = na.omit)
  ordinal <- as.numeric(model.response(mf))
  ellipsis.arg <- list(...)
  
  #print(head(data))
  
  if (is.null(ellipsis.arg$m)) {
    ellipsis.arg$m <- length(levels(factor(ordinal, 
                                           ordered = TRUE)))
  }
  
  if (is.null(ellipsis.arg$maxiter)) {
    ellipsis.arg$maxiter <- 500
  }
  if (is.null(ellipsis.arg$toler)) {
    ellipsis.arg$toler <- 1e-04
  }
  if (is.null(ellipsis.arg$iterc)) {
    ellipsis.arg$iterc <- 5
  }
  
  if (is.null(ellipsis.arg$invgen)) {
    ellipsis.arg$invgen<-TRUE
  }
  
  if (is.null(ellipsis.arg$starting)) {
    ellipsis.arg$invgen<-TRUE
  }
  
  covpai<-model.matrix(Formula,data=mf,rhs=1)
  covcsi<-model.matrix(Formula,data=mf,rhs=2)
  
  
  ##########
  if (ncol(covpai)==0){
    Y<-NULL
    if (ncol(covcsi)==0){
      W<-NULL
      dd<-data.frame('ordinal'=ordinal)
    } else {
      nomi_xi<-colnames(covcsi)[-1]
      W<-as.matrix(covcsi[,-1])
      if (is.null(nomi_xi)){
        colnames(W)<-paste("W",1:NCOL(W),sep="_")
      } else {
        colnames(W)<-nomi_xi
      }
      dd<-data.frame('ordinal'=ordinal,W)
      colnames(dd)<-c("ordinal",colnames(W))
      
    }
    
  } else {
    Y<-as.matrix(covpai[,-1])
    nomi_pai<-colnames(covpai)[-1]
    if (is.null(nomi_pai)){
      colnames(Y)<-paste("Y",1:NCOL(Y),sep="_")
    } else {
      colnames(Y)<-nomi_pai
    }
    
    if (ncol(covcsi)==0){
      W<-NULL
      dd<-data.frame('ordinal'=ordinal,Y)
      colnames(dd)<-c("ordinal",colnames(Y))
      
    } else {
      W<-as.matrix(covcsi[,-1])
      nomi_xi<-colnames(covcsi)[-1]
      
      if (is.null(nomi_xi)){
        colnames(W)<-paste("W",1:NCOL(W),sep="_")
      } else {
        colnames(W)<-nomi_xi
      }
      
      YW<-cbind(Y,W); nomiYW<-colnames(YW)
      qualicov<-unique(nomiYW)
      dd<-data.frame('ordinal'=ordinal,YW[,match(qualicov,nomiYW)])
      colnames(dd)<-c("ordinal",qualicov)
    }
    
  }
  
  ellipsis.arg$data<-dd
  
  lista<-ellipsis.arg
  
  m<-lista$m
  maxiter<-lista$maxiter
  toler<-lista$toler
  iterc<-lista$iterc
  invgen<-lista$invgen
  starting<-lista$starting
  
  if(is.null(Y) & is.null(W)) {
    if(m <= 3) stop("Number of ordered categories should be at least 4")
    mod<-fastcub00(m,ordinal,starting,maxiter,toler,iterc,invgen)
    
  }
  
  else{
    if(!is.null(Y) & is.null(W)) {
      Y<-as.matrix(Y)
      mod<-fastcubp0(m,ordinal,Y,starting,maxiter,toler,iterc,invgen)
      
    }
    else{
      if(is.null(Y) & !is.null(W)) {
        W<-as.matrix(W)
        mod<-fastcub0q(m,ordinal,W,starting,maxiter,toler,iterc,invgen)
        
      }
      else{
        if(!is.null(Y) & !is.null(W)) {
          Y<-as.matrix(Y)
          W<-as.matrix(W)
          mod<-fastcubpq(m,ordinal,Y,W,starting,maxiter,toler,iterc,invgen)
          
        }
        else cat("Wrong variables specification")
      }
    }                            
  }
  
  
  
  stime<-mod$estimates
  durata<-mod$time
  loglik<-as.numeric(mod$loglik)
  niter<-mod$niter
  varmat<-mod$vmatLouis
  BIC<-as.numeric(mod$BIC)
  time<-mod$durata
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,
                'BIC'=BIC,'ellipsis'=ellipsis.arg,
                'formula'=Formula,'call'=call)
  
  attr(results,"hidden")<-c("estimates","ordinal","loglik","varmat","BIC","ellipsis","family")
  class(results)<-"fastCUB"
  return(results)
  
}





