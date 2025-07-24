#' @title S3 method: summary for class "fastCUB"
#' @description S3 method summary for objects of class \code{\link{fastCUB}}. 
#' @aliases summary.fastCUB
#' @method summary fastCUB
#' @param object An object of class \code{\link{fastCUB}}
#' @param correlation Logical: should the estimated correlation matrix be returned? Default is FALSE
#' @param ...  Other arguments
#' @export 
#' @return Extended summary results of the fitting procedure, including parameter estimates, their standard errors and
#' Wald statistics, maximized log-likelihood compared with that of the saturated model and of a Uniform sample.
#' AIC, BIC and ICOMP indeces are also displayed for model selection. Execution time and number of exectued iterations 
#' for the fitting procedure are aslo returned.
#' @import methods
#' @rdname summary.fastCUB
#' @keywords package



################################################################
###############################################################


#summary <- function(object,correlation=FALSE,...) UseMethod("summary", object)

# digits=options()$digits 

summary.fastCUB <- function(object, correlation=FALSE, ...){
  
  flagcov<-0
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  m<-ellipsis$m
  n<-length(object$ordinal)
  
  output<-list()
  
  stime<-object$estimates
  ordinal<-object$ordinal
  
  freq<-tabulate(ordinal,nbins=m)
  varm<- vcov(object)    #as.matrix(object$varmat);
  np<-length(stime)
  if (isTRUE(varm==matrix(NA,nrow=np,ncol=np))==TRUE){
    trvarmat<-output$ICOMP<-NA
    output$errstd<-output$wald<-output$pval<-rep(NA,np)
  } else {
    trvarmat<-sum(diag(varm))
    output$loglik<- as.numeric(logLik(object))
    output$ICOMP<- -2*output$loglik + np*log(trvarmat/np) - log(det(varm))
    output$errstd<-sqrt(diag(varm));
    output$wald<-stime/output$errstd;
    output$pval<-round(2*(1-pnorm(abs(output$wald))),digits)
  }
  
  output$loglik<-logLik(object)
  output$AIC<- -2*logLik(object)+2*(np)
  output$BIC<- -2*logLik(object)+log(n)*(np)
  output$ellipsis<-ellipsis
  output$llunif<- -n*log(m);
  
  nonzero<-which(freq!=0)
  output$logsat <- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  output$devian<-2*(output$logsat-object$loglik)
  output$object<-object
  output$n<-n
  
  output$cormat<-NULL
  if (correlation==TRUE){
    output$cormat<- cormat(object)
  }
  StdErr<-output$errstd
  Wald<-output$wald
  matout<-cbind(stime,StdErr,Wald)
  
  colnames(matout)<-c("Estimates","StdErr","Wald")
  rownames(matout)<-parnames.fastCUB(object)
  output$results<-matout
  
  class(output)<-"summary.fastCUB"
  
  
  
  print.summary.fastCUB <- function(x,...){
    
    if(!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl, control = NULL)
    }
    
    ellipsis<-x$ellipsis
    object<-x$object
    maxiter<-object$ellipsis$maxiter
    family<-object$family
    niter<-object$niter
    m<-object$ellipsis$m
    stime<-object$estimates
    
    loglik<-   x$loglik
    aic<-x$AIC
    bic<-x$BIC
    
    modello<-object$formula
    data<-ellipsis$data
    
    mf<-model.frame(modello,data=data,na.action=na.omit)
    
    n<-x$n
    cat("=======================================================================","\n")
    cat("=====>>>", family," model    <<<=====   ML-estimates via fast EM algorithm  ","\n")
    cat("=======================================================================","\n")
    cat(" m=", m," Sample size: n=",n," Iterations=", niter," Maxiter=",maxiter,"\n")
    cat("=======================================================================","\n")
    cat(" Loglik =", loglik," AIC = ",aic," BIC =",bic,"\n")
    cat("=======================================================================","\n")
    
    StdErr<-x$errstd
    Wald<-x$wald
    
    
    data<-object$data
    
    listanomi<-parnames.fastCUB(object)
    
    covpai<-model.matrix(modello,data=mf,rhs=1)
    covcsi<-model.matrix(modello,data=mf,rhs=2)
    
    
    if (ncol(covpai)!=0 | ncol(covcsi)!=0 ){
      flagcov<-1
    }
    
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
    
    
    if ( !is.null(Y) & !is.null(W)){
      Y<-as.matrix(Y); W<-as.matrix(W); 
      p<-NCOL(Y);
      q<-NCOL(W); 
      
      mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
      colnames(mat1)<-c("Estimates","StdErr","Wald")
      rownames(mat1)<-listanomi[1:(p+1)]
      x$uncertainty<-mat1
      
      mat2<-cbind(stime[(p+2):(p+q+2)],StdErr[(p+2):(p+q+2)],Wald[(p+2):(p+q+2)])
      x$feeling<-mat2
      
      
      cat("Uncertainty                                           ", "\n")
      
      print(mat1,digits=digits)
      cat("=======================================================================","\n")
      cat("Feeling                                           ", "\n")
      colnames(mat2)<-c("Estimates","StdErr","Wald")
      rownames(mat2)<-listanomi[(p+2):(p+q+2)]
      
      print(mat2,digits=digits)
      cat("=======================================================================","\n")
      
      
    }  else if (is.null(Y) & !is.null(W)){
      W<-as.matrix(W);
      
      q<-NCOL(W); 
      
      mat1<-cbind(stime[1],StdErr[1],Wald[1])
      colnames(mat1)<-c("Estimates","StdErr","Wald")
      rownames(mat1)<-listanomi[1]
      
      mat2<-cbind(stime[2:(q+2)],StdErr[2:(q+2)],Wald[2:(q+2)])
      colnames(mat2)<-c("Estimates","StdErr","Wald")
      rownames(mat2)<-listanomi[2:(q+2)]
      cat("Uncertainty                                           ", "\n")
      
      print(mat1,digits=digits)
      cat("=======================================================================","\n")
      
      cat("Feeling                                           ", "\n")
      print(mat2,digits=digits)
      x$uncertainty<-mat1
      x$feeling<-mat2
      cat("=======================================================================","\n")
      
    } else if ( !is.null(Y) & is.null(W)){
      Y<-as.matrix(Y);
      
      p<-NCOL(Y); 
      
      mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
      colnames(mat1)<-c("Estimates","StdErr","Wald")
      rownames(mat1)<-listanomi[1:(p+1)]
      
      mat2<-cbind(stime[(p+2)],StdErr[p+2],Wald[p+2])
      colnames(mat2)<-c("Estimates","StdErr","Wald")
      rownames(mat2)<-listanomi[p+2]
      cat("Uncertainty                                           ", "\n")
      
      print(mat1,digits=digits)
      cat("=======================================================================","\n")
      
      cat("Feeling                                           ", "\n")
      print(mat2,digits=digits)
      x$uncertainty<-mat1
      x$feeling<-mat2
      cat("=======================================================================","\n")
      
    } else if ( is.null(Y) & is.null(W)) {
      
      mat1<-cbind(stime[1],StdErr[1],Wald[1])
      colnames(mat1)<-c("Estimates","StdErr","Wald")
      
      mat2<-cbind(stime[2],StdErr[2],Wald[2])
      colnames(mat2)<-c("Estimates","StdErr","Wald")
      rownames(mat2)<-listanomi[2]
      x$uncertainty<-mat1
      x$feeling<-mat2
      
      cat("Uncertainty                                           ", "\n")
      
      print(mat1,digits=digits)
      cat("=======================================================================","\n")
      
      cat("Feeling                                           ", "\n")
      print(mat2,digits=digits)
      
      cat("=======================================================================","\n")
      
    }
    
    
  }
  
  print(output)
  invisible(output$results)
  
  
}
