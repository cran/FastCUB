#' @title S3 method  vcov() for class "fastCUB"
#' @description S3 method: vcov for objects of class \code{\link{fastCUB}}. 
#' @aliases vcov.fastCUB
#' @param object An object of class \code{\link{fastCUB}}
#' @param ...  Other arguments
#' @method vcov fastCUB
#' @export 
#' @return Variance-covariance matrix of the final ML estimates for parameters of the fitted CUB model.
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
  
  ellipsis<-object$ellipsis
  family<-object$family
  varcov<-as.matrix(object$varmat)
  listanomi<-parnames.fastCUB(object)
  
  if (NROW(varcov)>1){
    rownames(varcov)<-colnames(varcov)<-listanomi
    #dimnames(varcov)<-list(listanomi,listanomi)
  } else {
    rownames(varcov)<-listanomi
    colnames(varcov)<-"Squared Standard Error"
    #dimnames(varcov)<-list(listanomi,"Squared Standard Error")
  }
  return(round(varcov,digits=digits))
}
