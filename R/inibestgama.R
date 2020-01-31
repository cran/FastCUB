#' @title Preliminary parameter estimates of a CUB model with covariates for feeling
#' @description Compute preliminary parameter estimates for the feeling component of a CUB model 
#' fitted to ordinal responses
#' These estimates are set as initial values for parameters to start the E-M algorithm.
#' @aliases inibestgama
#' @usage inibestgama(m,ordinal,W)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses 
#' @param W Matrix of selected covariates for explaining the feeling component
#' @export inibestgama
#' @return A vector of length equal to NCOL(W)+1, whose entries are the preliminary estimates
#'  of the parameters for the feeling component, including an intercept term as first entry.
#' @references Iannario M. (2008). Selecting feeling covariates in rating surveys, 
#' \emph{Rivista di Statistica Applicata}, \bold{20}, 103--116 \cr
#' Iannario M. (2009). A comparison of preliminary estimators in a class of ordinal data models,
#' \emph{Statistica & Applicazioni}, \bold{VII}, 25--44 \cr
#' Iannario M. (2012).  Preliminary estimators for a mixture model of ordinal data, 
#' \emph{Advances in Data Analysis and Classification}, \bold{6}, 163--184
#' @seealso \code{\link{inibest}}
#' @keywords htest utilities
#' @examples
#' data(univer)
#' m<-7; ordinal<-univer$global; cov<-univer$gender
#' ini<-inibestgama(m,ordinal,W=cov)


inibestgama<-function(m,ordinal,W){
  
  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  W<-as.matrix(W)
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  
  WW<-cbind(1,W)                           
  ni<-log((m-ordinal+0.5)/(ordinal-0.5))
  rr<-  qr(t(WW)%*%WW)$rank
  
  if (rr==ncol(WW)){
    mat<-solve(t(WW)%*%WW)
    gama<-(mat)%*%(t(WW)%*%ni) 
  } else {
    gama<-as.matrix(rep(0.1,ncol(WW)))
  }
  
  q<-NCOL(W)
  listanomi<-paste("gamma",0:q,sep="_")
  dimnames(gama)<-list(listanomi,"")
  return(gama)
}
