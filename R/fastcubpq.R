#' @title Main function for CUB models with covariates for both the uncertainty and the feeling components
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining both the
#'  feeling and the uncertainty components by means of logistic transform.
#' @aliases fastcubpq
#' @usage fastcubpq(m,ordinal,Y,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param W Matrix of selected covariates for explaining the feeling component
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm
#' @param toler Fixed error tolerance for final estimates
#' @param iterc Iteration from which the acceleration strategy starts
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @param verbose Logical: should messages about  acceleration steps be printed out? (Default is FALSE)
#' @return An object of the class "fastCUB"
#' @import stats
#' @seealso \code{\link{loglikcubpq}}, \code{\link{inibestgama}}
#' @keywords internal

########################################

########################################
########################################
########################################
fastcubpq<-function(m,ordinal,Y,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE){
  
  tt0<-proc.time()
  n<-length(ordinal)
  Y<-as.matrix(Y)
  W<-as.matrix(W)
  p<-NCOL(Y)
  q<-NCOL(W)
  aver<-mean(ordinal)
  
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  
  YY<-cbind(1,Y);     WW<-cbind(1,W);          
  #################################################################################
  freq<-tabulate(ordinal,nbins=m)
  
  if (!is.null(starting)){
    betjj<-starting[1:(p+1)]; gamajj<-starting[(p+2):(p+q+2)]
  } else {
    inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; bet0<-log(pai/(1-pai)); betjj<-c(bet0,rep(0.1,p));
    gamajj<-as.numeric(inibestgama(m,ordinal,W) )
  }
  
  #################################################################################
  loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,q) ***************************
  # ********************************************************************
  nniter<-1
  lb1<-rep(-5,p+1); ub1<-rep(5,p+1)
  lb2<-rep(-5,q+1); ub2<-rep(5,q+1)
  iterc2<-iterc
  
  while(nniter<=maxiter){
    
    thetaold<-c(betjj,gamajj)
    
    loglikold<-loglikjj
   
    vettn<-bitgama(m,ordinal,W,gamajj)   
    
    
    aai<- -1+1/(logis(Y,betjj))  
    ttau<-1/(1+aai/(m*vettn))  
    ####################  maximize w.r.t. bet and gama    ############
    esterno10<-cbind(ttau,YY)
    
    esterno01<-cbind(ttau,ordinal,WW)
    bet<-betjj;  gama<-gamajj;
    
    betoptim<-optim(bet,effe10,esterno10=esterno10)
    gamaoptim<-optim(gama,effe01,esterno01=esterno01,m=m)
    
   
    # ################################################################         
    betjj<-betoptim$par
    gamajj<-gamaoptim$par
    
    thetanew<-c(betjj,gamajj)
    param<-thetanew
    loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
    
    
    if (nniter>= iterc2){
      
      csijj<-logis(W,gamajj)
      ai<- ordinal - 1 -  (m-1)*(1-csijj)
      vettn<-as.numeric(bitgama(m,ordinal,W,gamajj)   )
      aai<- -1+1/(logis(Y,betjj))   
      ttau<-1/(1+aai/(m*vettn))        
      
      dc<-decomp(ttau,ordinal,m,param,ai,Y=Y,W=W)
      
      M1<-dc$vcScorec
      M2<-dc$vcScore
      InfC<-dc$Ic
      
      D<-matrix(0,nrow=length(param),ncol=length(param))
      diag(D)<-diag(InfC)
      D2<-InfC - D 
      
      H<- - M1 + M2 + D2
      
      listE<-list()
      
      for (l in 1:nrow(H)){
        listE[[l]]<-matrix(0,nrow=length(param),ncol=length(param))
        listE[[l]][l,]<-H[l,]
      }
      
      if (invgen==TRUE){
        Iinv<- invmatgen(D,H,listE);
      } else {
        Iinv<-try(solve(D+H),silent=TRUE)
      }
      
      if (any(class(Iinv)=="try-error")){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if (any(is.infinite(Iinv)) | any(is.nan(Iinv))){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if (all(class(Iinv)!="try-error") && !any(is.na(Iinv))){
        #
        if (verbose==TRUE) {print("accelerating")}
        
        loglikold2<-loglikjj
        
        inv<-Iinv%*%InfC
        dk<- thetanew-thetaold
        
        thetastar<- thetanew + inv%*%dk
        
        betjj2<-thetastar[1:(p+1)]; gamajj2<-thetastar[(p+2):(p+q+2)]
        loglikjj2<-loglikcubpq(m,ordinal,Y,W,betjj2,gamajj2)
        
        if ( is.nan(loglikjj2)) {
          Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
          if (verbose==TRUE) {cat("not accelerating at iter",nniter,"\n")}
        
          nniter <- nniter + 1
          next
        } else if (loglikjj2 < loglikold2){
          Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
          if (verbose==TRUE) {cat("not accelerating at iter",nniter,"\n")}
        
          nniter <- nniter + 1
          next
        } else {
          betjj<-betjj2; gamajj<-gamajj2; 
          loglikold<-loglikold2
          loglikjj<-loglikjj2
        }
        
        

        
      }  else {
        if (verbose==TRUE) {cat("not accelerating at iter",nniter,"\n")}
        #Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
        iterc2<-iterc2 + 10
        nniter <- nniter + 1
        next
      }
      
      
      
    }  
    
    
    
    testll<-abs(loglikjj-loglikold)
    
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  
  
  if (nniter < iterc ){
    
    param<-thetanew
    csijj<-logis(W,gamajj)
    ai<- ordinal - 1 -  (m-1)*(1-csijj)
    vettn<-as.numeric(bitgama(m,ordinal,W,gamajj)   )
    aai<- -1+1/(logis(Y,betjj))   
    ttau<-1/(1+aai/(m*vettn))        
    
    
    dc<-decomp(ttau,ordinal,m,param,ai,Y=Y,W=W)
    
    M1<-dc$vcScorec
    M2<-dc$vcScore
    InfC<-dc$Ic
    
    D<-matrix(0,nrow=length(param),ncol=length(param))
    diag(D)<-diag(InfC)
    D2<-InfC - D 
    
    H<- - M1 + M2 + D2
    
    listE<-list()
    
    for (l in 1:nrow(H)){
      listE[[l]]<-matrix(0,nrow=length(param),ncol=length(param))
      listE[[l]][l,]<-H[l,]
    }
    
    if (invgen==TRUE){
      Iinv<- invmatgen(D,H,listE)
    } else {
      Iinv<- try(solve(D+H),silent=TRUE)
    }
    if (any(class(Iinv)=="try-error")){
      Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
    }
    
  }
  
  durata<-proc.time()-tt0;durata<-durata[1];
 
  if (any(is.infinite(Iinv)) || any(is.nan(Iinv))){
    Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
  }
  
  if(any(is.na(Iinv))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-vmatLouis<-matrix(NA,nrow=nrow(Iinv),ncol=nrow(Iinv))
    se<-rep(NA,nrow(Iinv))
  } else {
    if(det(Iinv)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-vmatLouis<-matrix(NA,nrow=nrow(Iinv),ncol=nrow(Iinv))
      se<-rep(NA,nrow(Iinv))
    } else {
      varmat<-vmatLouis<-Iinv
      se<-sqrt(diag(varmat))
    }
  }
  
  
  
  
  bet<-betjj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUBpq<- -2*loglik+2*(p+q+2)
  BICCUBpq<- -2*loglik+log(n)*(p+q+2)
  
  stime<-c(bet,gama)
  
  
  results<-list('estimates'=stime,'se'=se,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBpq,'vmatLouis'=vmatLouis)
  
  return(results)
}
