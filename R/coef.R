#' @title S3 Method: coef for class "fastCUB"
#' @description S3 method: coef for objects of class \code{\link{fastCUB}}. 
#' @aliases coef.fastCUB
#' @method coef fastCUB
#' @param object An object of class \code{\link{fastCUB}}
#' @param ...  Other arguments
#' @import methods
#' @return ML estimates of parameters of the fitted CUB model.
#' @details Returns estimated values of coefficients of the fitted model 
#' @export 
#' @rdname coef.fastCUB
#' @keywords package
#' @seealso \code{\link{fastCUB}}, \code{\link{summary}}
#coef<- function(object,...) UseMethod("coef", object)


coef.fastCUB<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  output<-list()
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  listanomi<-parnames(object)
  mat<-round(as.matrix(object$estimates),digits=digits)
  dimnames(mat)<-list(listanomi,"")
  mat
  
}




