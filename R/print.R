#' @title S3 method: print for class "fastCUB"
#' @description S3 method print for objects of class \code{\link{fastCUB}}. 
#' @aliases print.fastCUB
#' @method print fastCUB
#' @param x An object of class \code{\link{fastCUB}}
#' @param ...  Other arguments
#' @export 
#' @return Brief summary results of the fitting procedure, including parameter estimates, their standard errors and 
#' the executed call.
#' @import methods
#' @rdname print.fastCUB
#' @keywords package



print.fastCUB<-function(x,...){
  
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  
  sterr<-as.numeric(round(x$se,digits=digits))
  
  mat<-cbind(round(x$estimates,digits=digits),sterr)
  rownames(mat)<-parnames(x)
  colnames(mat)<-c("Estimates","Standard Errors")
  
  object<-x$object
  stime<-object$estimates
  cat("","\n")
  print(mat)

  cat("","\n")
  cat("Maximized Log-Likelihood:",logLik(x,digits=digits),"\n")
  
  invisible(x)
}














