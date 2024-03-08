#' @title Generic function for coefficient names 
#' @description Generic function for names of parameter estimates of object of class "fastCUB". 
#' @aliases parnames.fastCUB
#' @param object An object of class "fastCUB"
#' @import methods
#' @return Parameter names
#' @keywords internal
#' @seealso \code{\link{summary}}

parnames <- function(object) UseMethod("parnames", object)




parnames.fastCUB<-function(object){
  
  effe<-object$formula
  #   EFFE<-modello$Formula
  data<-object$ellipsis$data
  mf<-model.frame(effe,data=data,na.action=na.omit)
  
  covpai<-model.matrix(effe,data=mf,rhs=1)
  covcsi<-model.matrix(effe,data=mf,rhs=2)

  if (ncol(covpai)==0){
    Y<-NULL
  } else {
    
    if (NCOL(covpai)==2){
      Y<-as.matrix(covpai[,-1])
      colnames(Y)<-colnames(covpai)[2]
    } else {
      Y<-covpai[,-1]
    }
  }
  if (ncol(covcsi)==0){
    W<-NULL
  } else {
    if (NCOL(covcsi)==2){
      W<-as.matrix(covcsi[,-1])
      colnames(W)<-colnames(covcsi)[2]
    } else {
      W<-covcsi[,-1]
    }
  }

  
 
  ellipsis<-object$ellipsis
  listanomi<-c()

    if (is.null(Y) & is.null(W)){
      listanomi<-c("pai","csi")
    } 
     if (!is.null(Y) & is.null(W)){
      betacoef<-c()
      npar<-length(object$estimates)
      if (!is.null(colnames(Y))){
        betacoef<-c("constant",colnames(Y))
      } else {
      
       for (j in 1:(npar-1)){
        betacoef[j]<-paste("beta",j-1,sep="_")
       }
      }
       listanomi<-c(betacoef,"csi")
     } 
    
    if (is.null(Y) & !is.null(W)){
      gamacoef<-c()
      npar<-length(object$estimates)
     
      if (!is.null(colnames(W))){
        gamacoef<-c("constant",colnames(W))
      } else {
        for (j in 1:(npar-1)){
          gamacoef[j]<-paste("gamma",j-1,sep="_")
        }
      }
      
      listanomi<-c("pai",gamacoef)
    }
    
    if (!is.null(Y) & !is.null(W)) {
      betacoef<-gamacoef<-c()
      Y<-as.matrix(Y); W<-as.matrix(W)
      ny<-NCOL(Y); nw<-NCOL(W);
      
      if (is.null(colnames(Y))){
        for (j in 1:(ny+1)){
          betacoef[j]<-paste("beta",j-1,sep="_")
        }
        
        } else {
          betacoef<-c("constant",colnames(Y))
        }
        
      
      
      if (is.null(colnames(W))){
        for (j in 1:(nw+1)){
          gamacoef[j]<-paste("gamma",j-1,sep="_")
        }
        
      } else {
        gamacoef<-c("constant",colnames(W))
        
      }
     
      listanomi<-c(betacoef,gamacoef)
    }
    
  
  return(listanomi)
}

