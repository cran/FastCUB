#' @title S3 BIC method for class "fastCUB"
#' @description S3 BIC method for objects of class \code{\link{fastCUB}}. 
#' @aliases BIC.fastCUB
#' @method BIC fastCUB
#' @param object An object of class "fastCUB"
#' @param ...  Other arguments
#' @export 
#' @return BIC index for the fitted model.
#' @rdname BIC.fastCUB
#' @import methods
#' @seealso \code{\link{logLik}}, \code{\link{fastCUB}}
#' @keywords package

#BIC<- function(object,...) UseMethod("BIC", object)

BIC.fastCUB<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  bic<-object$BIC
  round(bic,digits=digits)
}

