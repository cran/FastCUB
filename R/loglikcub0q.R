#' @title Log-likelihood function of a CUB model with covariates for the feeling component
#' @description Compute the log-likelihood function of a CUB model fitting ordinal data, with \eqn{q} 
#' covariates for explaining the feeling component.
#' @aliases loglikcub0q
#' @usage loglikcub0q(m, ordinal, W, beta0, gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates for explaining the feeling component
#' @param beta0 Logit transform of uncertainty parameter
#' @param gama Vector of parameters for the feeling component, with length NCOL(W) + 1 to account for 
#' an intercept term (first entry of gama)
#' @keywords internal
#' 


loglikcub0q <-
  function(m,ordinal,W,beta0,gama){
    pai<-1/(1+exp(-beta0))
    W<-as.matrix(W)
    probn<-probcub0q(m,factor(ordinal,ordered=TRUE),W,pai,gama)
    sum(log(probn))
  }

