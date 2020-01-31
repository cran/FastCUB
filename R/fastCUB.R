#' @title Main function for fast estimation CUB models
#' @description Main function to estimate and validate a CUB model for explaining uncertainty
#' and feeling for given ratings, with or without covariates, on the basis of Louis' identity for the information matrix and the derived accelerated estimation.
#' @aliases fastCUB
#' @usage fastCUB(Formula, data, mix=FALSE, tolmix=1e+2,fmix= NULL,...)
#' @param Formula Object of class Formula with two right-hand side: the first for uncertainty covariates, the second for feeling covariates.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param mix Logical: should a first preliminary standard EM be run at toler equal to tolmix? (default is FALSE).
#' @param tolmix Convergence tolerance for first preliminary EM (if mix=TRUE).
#' @param fmix Fraction of iteration needed for preliminary EM (if mix=TRUE). Default is null.
#' @param ... Additional arguments to be passed for the specification of the model and the acceleration steps.
#' @return An object of the class "fastCUB": returns a list containing the following results:
#' \item{estimates}{Maximum likelihood estimates of model parameters}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' \item{parnames}{Names for model parameters}
#' @import Formula methods
#' @keywords package
#' @export fastCUB
#' @details This is the main function for CUB models, which calls for the corresponding functions whenever
#' covariates are specified. It performs maximum likelihood estimation via the E-M algorithm
#' for CUB models and extensions based on the Louis'identity for the observed information matrix.
#' @seealso \code{\link{probcub00}}, \code{\link{probcubp0}}, \code{\link{probcub0q}}, \code{\link{probcubpq}}
#' @examples
#' \donttest{
#' library(FastCUB)
#' data(univer)
#' ordinal<-univer$global
#' m<-7
#' effe<-with(univer, Formula(global~0|gender+freqserv+age +changefa))
#' cub0q<-fastCUB(effe,data=univer,m=7, maxiter=100,toler=1e-8,mix=TRUE,verbose=FALSE)
#' summary(cub0q)
#' ## Fast EM for  CUB model with covariates only for uncertainty
#' effe<-with(univer, Formula(global~gender+freqserv+age +changefa|0))
#' cubp0<-fastCUB(effe,data=univer,m=7, maxiter=100,toler=1e-8,iterc=5,verbose=TRUE)
#' ## Fast EM for  CUB model with covariates for both feeling and uncertainty
#' effe<-with(univer, Formula(global~gender+freqserv+age +changefa|gender+freqserv+age +changefa))
#' cubpq<-fastCUB(effe,data=univer,m=7, maxiter=100,toler=1e-8,iterc=5)
#' summary(cubpq)
#' BIC(cubpq)
#' }
#'

fastCUB<-function(Formula,data,mix=FALSE, tolmix=1e+2,fmix= NULL,...){


  call <- match.call()
  if (missing(data)){
    data<-environment(Formula)
  }

  mf <- model.frame(Formula, data, na.action = na.omit)

  ordinal <- as.numeric(model.response(mf))
  ellipsis.arg <- list(...)
  
  if (is.null(ellipsis.arg[["m"]])) {
    ellipsis.arg[["m"]] <- length(levels(factor(ordinal, ordered = TRUE)))
  }

  if (is.null(ellipsis.arg[["maxiter"]])) {
    ellipsis.arg[["maxiter"]] <- 500
  }
  if (is.null(ellipsis.arg$toler)) {
    ellipsis.arg$toler <- 1e-04
  }
  if (is.null(ellipsis.arg$invgen)) {
    ellipsis.arg$invgen <- TRUE
  }
  if (is.null(ellipsis.arg$iterc)) {
    ellipsis.arg$iterc <- 3
  }
  #
  if (is.null(ellipsis.arg$verbose)) {
    ellipsis.arg$verbose <- FALSE
  }
  covpai<-model.matrix(Formula,data=mf,rhs=1)
  covcsi<-model.matrix(Formula,data=mf,rhs=2)


  if (NCOL(covpai)==0){
    Y<-NULL
  } else {
    Y<-as.matrix(covpai[,-1])
    colnames(Y)<-colnames(covpai)[-1]
  }
  if (NCOL(covcsi)==0){
    W<-NULL
  } else {
    W<-as.matrix(covcsi[,-1])
    colnames(W)<-colnames(covcsi)[-1]
  }


  if (!is.null(Y) & !is.null(W)){
    
    effecub<-Formula(ordinal ~ Y|W|0)
  } else {
    if (is.null(Y)){
     
      effecub<-Formula(ordinal ~ 0|W|0)
    } else {
    
      effecub<-Formula(ordinal ~ Y|0|0)
    }
  }
  
  

  lista<-ellipsis.arg;

  m<-lista[['m']]
  maxiter<-lista[['maxiter']]
  toler<-lista[['toler']]
  iterc<-lista[['iterc']]
  invgen<-lista[['invgen']]
  verbose<-lista[['verbose']]
  starting<-lista[['starting']]

  
  if (mix == TRUE) {
    requireNamespace("CUB", "GEM")
    
    cub <- try(GEM(effecub, family = "cub", toler = tolmix),silent=TRUE)
    if (all(class(cub)!="try-error")){
      iterc <- ifelse(!is.null(fmix), round(fmix*cub$niter), 1)
      starting <- cub$estimates
      nitercub<-cub$niter; timecub<-cub$time
    } else {
      starting = NULL
      iterc<-10
      nitercub<-0; timecub<-0
    }
  }  else {
    starting = NULL
  }
  
  
  
  if(is.null(Y) & is.null(W)) {
    if(m <= 3) stop("Number of ordered categories should be at least 4")
    mod<-try(fastcub00(m,ordinal,starting,maxiter,toler,iterc,invgen,verbose),silent=TRUE)
    listanomi<-c("pai","csi")
    
  }

  else{
    if(!is.null(Y) & is.null(W)) {
      Y<-as.matrix(Y)
      mod<-try(fastcubp0(m,ordinal,Y,starting,maxiter,toler,iterc,invgen,verbose),silent=TRUE)
      listanomi<-c("intercept",colnames(Y),"csi")
    }
    else{
      if(is.null(Y) & !is.null(W)) {
        W<-as.matrix(W)
        mod<-try(fastcub0q(m,ordinal,W,starting,maxiter,toler,iterc,invgen,verbose),silent=TRUE)
       listanomi<-c("pai","intercept",colnames(W))

      }
      else{
        if(!is.null(Y) & !is.null(W)) {
          Y<-as.matrix(Y)
          W<-as.matrix(W)
          listanomi<-c("intercept",colnames(Y),"intercept",colnames(W))
          mod<-try(fastcubpq(m,ordinal,Y,W,starting,maxiter,toler,iterc,invgen,verbose),silent=TRUE)

        }
        else warning("Wrong variables specification")
      }
    }
  }

  if (class(mod)=="try-error"){
    warning("Unable to estimate model: try to change tolerance and estimation settings","\n")
    stime<-NA
    durata<-NA
    loglik<-NA
    niter<-NA
    varmat<-NA
    se<-NA
    BIC<-Inf
    time<-NA
    envc<-NA
  } else {
    stime<-mod$estimates
    durata<-mod$time + ifelse(mix==TRUE,timecub,0)
    loglik<-as.numeric(mod$loglik)
    niter<-mod$niter + ifelse(mix==TRUE,nitercub,0)
    varmat<-mod$vmatLouis
    se<-mod$se
    BIC<-as.numeric(mod$BIC)
    time<-mod$durata
    
    envc<-new.env()
    envc$data<-data
    envc$formula<-Formula
    envc$covpai<- covpai
    envc$covcsi<-covcsi
  }



  results<-list('estimates'=stime,'se'=se,'env'=envc,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,'parnames'=listanomi,
                'BIC'=BIC,'ellipsis'=ellipsis.arg,'covpai'=covpai,'covcsi'=covcsi,'formula'=Formula,'call'=call,'data'=data)

  attr(results,"hidden")<-c("estimates","ordinal","loglik","varmat","BIC","ellipsis")
  class(results)<-"fastCUB"
  return(results)

}











