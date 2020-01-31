#' @title Log-likelihood function of a CUB model without covariates
#' @description Compute the log-likelihood function of a CUB model without covariates for a given 
#' absolute frequency distribution.
#' @aliases loglikcub00
#' @usage loglikcub00(m, freq, beta0, gama0)
#' @param m Number of ordinal categories
#' @param freq Vector of the absolute frequency distribution
#' @param beta0 Logit transform of uncertainty parameter
#' @param gama0 Logit transform of fFeeling parameter
#' @keywords internal
#' 



loglikcub00 <-
  function(m,freq,beta0,gama0){
    pai<-1/(1+exp(-beta0)); csi<-1/(1+exp(-gama0))
    
    t(freq)%*%log(probcub00(m,pai,csi))
    
  }