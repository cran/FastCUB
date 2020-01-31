#' @title S3 method  vcov() for class "fastCUB"
#' @description S3 method: vcov for objects of class \code{\link{fastCUB}}. 
#' @aliases vcov.fastCUB
#' @param object An object of class \code{\link{fastCUB}}
#' @param ...  Other arguments
#' @method vcov fastCUB
#' @export 
#' @return Variance-covariance matrix of the final ML estimates for parameters of the fitted CUB model, 
#' according to the logit transform also when covariates are not included for uncertainty or feeling parameters. 
#' It is computed on the basis of Louis' identity within the EM algorithm.
#' @import methods
#' @seealso \code{\link{fastCUB}}
#' @keywords package
#' @rdname vcov.fastCUB



#vcov <- function(object,...) UseMethod("vcov", object)


vcov.fastCUB<-function(object, ...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
   varcov<-object$varmat
   listanomi<-parnames(object)
   if (NCOL(object$covpai)==0){
     listanomi[1]<-paste("beta","0",sep="_")
   } 
     np<-length(object$estimates)
     if (NCOL(object$covcsi)==0){
       listanomi[np]<-paste("gamma","0",sep="_")
     }
   
   
  
   rownames(varcov)<-listanomi
   colnames(varcov)<-listanomi
  
  return(round(varcov,digits=digits))
}
