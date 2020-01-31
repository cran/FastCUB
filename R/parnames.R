#' @title Generic function for coefficient names
#' @description Generic function for names of parameter estimates of object of class "fastCUB".
#' @aliases parnames
#' @param object An object of class "fastCUB"
#' @keywords internal
#' @return Parameter names
#' @seealso \code{\link{summary}}

#parnames <- function(object) UseMethod("parnames", object)


parnames<-function(object){


  covpai<-  object$covpai 
  covcsi<- object$covcsi 


  if (NCOL(covpai)==0){
    Y<-NULL
  } else {

    if (NCOL(covpai)==2){
      Y<-as.matrix(covpai[,-1])
    
    } else {
      Y<-covpai[,-1]
    }
  }
  if (NCOL(covcsi)==0){
    W<-NULL
  } else {
    if (NCOL(covcsi)==2){
      W<-as.matrix(covcsi[,-1])
    
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
      npar<-NCOL(Y) 
      if (!is.null(colnames(Y))){
        betacoef<-c("intercept",colnames(Y))
      } else {
        betacoef[1]<-"intercept"
       for (j in 2:(npar+1)){
        betacoef[j]<-paste("Y",j-1,sep="_")
       }
      }
       listanomi<-c(betacoef,"csi")
     }

    if (is.null(Y) & !is.null(W)){
      gamacoef<-c()
      npar<-NCOL(W) 

      if (!is.null(colnames(W))){
        gamacoef<-c("intercept",colnames(W))
      } else {
        gamacoef[1]<-"intercept"
        for (j in 2:(npar+1)){
          gamacoef[j]<-paste("W",j-1,sep="_")
        }
      }

      listanomi<-c("pai",gamacoef)
    }

    if (!is.null(Y) & !is.null(W)) {
      betacoef<-gamacoef<-c()
      Y<-as.matrix(Y); W<-as.matrix(W)
      ny<-NCOL(Y); nw<-NCOL(W);

      if (is.null(colnames(Y))){
        betacoef[1]<-"intercept"
        for (j in 2:(ny+1)){
          betacoef[j]<-paste("Y",j-1,sep="_")
        }

        } else {
          betacoef<-c("intercept",colnames(Y))
        }



      if (is.null(colnames(W))){
        gamacoef[1]<-"intercept"
        for (j in 2:(nw+1)){
          gamacoef[j]<-paste("W",j-1,sep="_")
        }

      } else {
        gamacoef<-c("intercept",colnames(W))

      }

      listanomi<-c(betacoef,gamacoef)
    }
  return(listanomi)
}

