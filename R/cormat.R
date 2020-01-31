#' @title Correlation matrix for estimated model
#' @description Compute parameter correlation matrix for estimated model as returned by an object
#'  of class "fastCUB". 
#' @aliases cormat
#' @usage cormat(object,digits=options()$digits)
#' @param object An object of class "fastCUB"
#' @param digits Number of significant digits to be printed. Default is \code{options()$digits}
#' @export 
#' @return Parameters correlation matrix for fitted fastCUB models. 
#' @keywords models
#' @seealso \code{fastCUB}, \code{vcov}


cormat<-function(object,digits=options()$digits){
  
  varmat<-vcov(object)
    np<-NROW(varmat)
    cormat<-diag(rep(1,np))

    if (np>1){
      if (isTRUE(varmat==matrix(NA,nrow=np,ncol=np))==TRUE){
        cormat<-matrix(NA,nrow=np,ncol=np)
        dimnames(cormat)<-list(parnames(object),parnames(object))

      } else {
        ddd<-diag(sqrt(1/diag(varmat)))
        cormat<-round(ddd%*%varmat%*%ddd,digits)
        dimnames(cormat)<-list(parnames(object),parnames(object))

      }

    } else {
      dimnames(cormat)<-list(parnames(object),"Correlation")
    }

    return(cormat)
  
  
}







# 
# cormat <- function(object) UseMethod("cormat", object)
# 
# 
# cormat.GEM<-function(object){
#   
#   if(!inherits(object, "GEM")) stop("not a \"GEM\" fit")
#   
#   
#   varmat<-vcov(object)
#   np<-NROW(varmat)
#   cormat<-diag(rep(1,np))
#   
#   if (np>1){
#     if (isTRUE(varmat==matrix(NA,nrow=np,ncol=np))==TRUE){
#       cormat<-matrix(NA,nrow=np,ncol=np)
#       dimnames(cormat)<-list(parnames(object),parnames(object))
#       
#     } else {
#       ddd<-diag(sqrt(1/diag(varmat)))
#       cormat<-round(ddd%*%varmat%*%ddd,5)  
#       dimnames(cormat)<-list(parnames(object),parnames(object))
#       
#     }
#     
#   } else {
#     dimnames(cormat)<-list(parnames(object),"Correlation")
#   }
#   
#   return(cormat)
# }
# 
