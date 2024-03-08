#' @title Normalized dissimilarity measure
#' @description Compute the normalized dissimilarity measure between observed
#'  relative frequencies and estimated (theoretical) probabilities of a discrete distribution.
#' @usage dissim(proba,probb)
#' @aliases dissim
#' @param proba Vector of observed relative frequencies
#' @param probb Vector of estimated (theoretical) probabilities
#' @return Numeric value of the dissimilarity index, assessing the distance to a perfect fit.
#' @keywords univar
#' @export dissim
#' @examples 
#' proba<-c(0.01,0.03,0.08,0.07,0.27,0.37,0.17)
#' probb<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' dissim(proba,probb)


dissim <-
function(proba,probb){
  
  if (length(proba)==length(probb)){
    
    return(0.5*sum(abs(proba-probb)))
  } else {
    cat("Error: input vectors should have the same length","\n")
  }
  

  
  }
