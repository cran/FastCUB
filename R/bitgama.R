#' @title Shifted Binomial distribution with covariates
#' @description Return the shifted Binomial probabilities of ordinal responses where the feeling component 
#' is explained by covariates via a logistic link.
#' @aliases bitgama
#' @usage bitgama(m,ordinal,W,gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of covariates for the feeling component
#' @param gama Vector of parameters for the feeling component, with length equal to 
#' NCOL(W)+1 to account for an intercept term (first entry of \code{gama})
#' @export bitgama
#' @return A vector of the same length as \code{ordinal}, where each entry is the shifted Binomial probability for
#'  the corresponding observation and feeling value.
#' @seealso  \code{\link{logis}}, \code{\link{probcub0q}}, \code{\link{probcubpq}} 
#' @keywords distribution
#' @import stats
#' @examples 
#' n<-100
#' m<-7
#' W<-sample(c(0,1),n,replace=TRUE)
#' gama<-c(0.2,-0.2)
#' csivett<-logis(W,gama)
#' ordinal<-rbinom(n,m-1,csivett)+1
#' pr<-bitgama(m,ordinal,W,gama)

bitgama <-function(m,ordinal,W,gama){
  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  W <- as.matrix(W)
  
  if (ncol(W)==1){
    W<-as.numeric(W)
  }

  ci<- 1/(logis(W,gama))-1
  kkk(m,ordinal)*exp((ordinal-1)*log(ci)-(m-1)*log(1+ci))
}
