#' @title logLik S3 Method for class "fastCUB"
#' @description S3 method: logLik() for objects of class "fastCUB". 
#' @aliases logLik.fastCUB
#' @method logLik fastCUB
#' @param object An object of class "fastCUB"
#' @param ...  Other arguments
#' @export 
#' @return Log-likelihood at the final ML estimates for parameters of the fitted fastCUB model.
#' @import methods
#' @seealso  \code{\link{fastCUB}}
#' @keywords package
#' @rdname logLik.fastCUB



logLik.fastCUB<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  return(round(object$loglik,digits=digits))    
  
}