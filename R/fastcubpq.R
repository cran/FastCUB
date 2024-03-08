#' @title Main function for CUB models with covariates for both the uncertainty and the feeling components
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining both the
#'  feeling and the uncertainty components by means of logistic transform.
#' @aliases fastcubpq
#' @usage fastcubpq(m,ordinal,Y,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param W Matrix of selected covariates for explaining the feeling component
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @return An object of the class "fastCUB"
#' @import stats
#' @seealso \code{\link{loglikcubpq}}, \code{\link{inibestgama}}
#' @keywords internal 

########################################

########################################
########################################
########################################
fastcubpq<-function(m,ordinal,Y,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE){
  
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
    gamajj<-inibestgama(m,ordinal,W) 
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
    #print(nniter)
    thetaold<-c(betjj,gamajj)
    
    loglikold<-loglikjj
    #cat("loglik old:", loglikold,"\n")
    vettn<-as.numeric(bitgama(m,ordinal,W,gamajj)   )
    aai<- -1+1/(logis(Y,betjj))   
    ttau<-1/(1+aai/(m*vettn))        
    ####################  maximize w.r.t. bet and gama    ############
    esterno10<-cbind(ttau,YY)
    esterno01<-cbind(ttau,ordinal,WW)
    bet<-betjj;  gama<-gamajj;
    betoptim<-optim(bet,effe10,esterno10=esterno10)
    gamaoptim<-optim(gama,effe01,esterno01=esterno01,m=m)
    
    # betoptim<-optim(bet,effe10,esterno10=esterno10,method="L-BFGS-B",lower=lb1,upper=ub1,hessian=TRUE)
    # gamaoptim<-optim(gama,effe01,esterno01=esterno01,m=m,method="L-BFGS-B",lower=lb2,upper=ub2,hessian=TRUE)
    # 
    # ################################################################         
    betjj<-betoptim$par
    gamajj<-gamaoptim$par
    
    thetanew<-c(betjj,gamajj)
    
    #ai<- (m-1)*(1-csii) - (ordinal-1)
    #print(nniter)
    if (nniter>= iterc2){
      
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
        Iinv<- invmatgen(D,H,listE);
      } else {
        Iinv<-solve(D+H)
      }
      
      # print("accelerating")
      #     inv<-Iinv%*%InfC
      #     dk<- thetanew-thetaold
      # 
      #     thetastar<- thetanew + inv%*%dk
      #   #  loglikold<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
      #     betjj<-thetastar[1:(p+1)]; gamajj<-thetastar[(p+2):(p+q+2)]
      
      if (!is(Iinv,"try-error") && !any(is.nan(Iinv))){
        
        if (rcond(Iinv)> 1e-8){
          print("accelerating")
          inv<-Iinv%*%InfC
          dk<- thetanew-thetaold
          
          thetastar<- thetanew + inv%*%dk
          # loglikold<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
          betjj<-thetastar[1:(p+1)]; gamajj<-thetastar[(p+2):(p+q+2)]
          
        } else {
          cat("not accelerating at iter",nniter,"\n")
          Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
          iterc2<-iterc2 + 10
          nniter <- nniter + 1
          next
        }
        
      }  else {
        cat("not accelerating at iter",nniter,"\n")
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
        iterc2<-iterc2 + 10
        nniter <- nniter + 1
        next
      }
      
      
      
    }  ## chiude if nniter >= iterc2
    
    
    loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
    
    testll<-abs(loglikjj-loglikold)
    #print(testll)
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
      if (is(Iinv,"try-error")){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
    } else {
      Iinv<- solve(D+H)
    }
    
  }
  
  durata<-proc.time()-tt0;durata<-durata[1];
  # cat("iterc2",iterc2,"\n")
  #  cat("nniter",nniter,"\n")
  
  if(any(is.na(Iinv))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-vmatLouis<-matrix(NA,nrow=nrow(Iinv),ncol=nrow(Iinv))
  } else {
    if(det(Iinv)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-vmatLouis<-matrix(NA,nrow=nrow(Iinv),ncol=nrow(Iinv))
    } else {
      varmat<-vmatLouis<-Iinv
    }
  }
  
  bet<-betjj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUBpq<- -2*loglik+2*(p+q+2)
  BICCUBpq<- -2*loglik+log(n)*(p+q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  #varmat<-varcovcubpq(m,ordinal,Y,W,bet,gama)
  #if(det(varmat)<=0) stop("Variance-covariance matrix NOT positive definite")
  #nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),paste("gamma",0:(length(gama)-1),sep="_"))
  stime<-c(bet,gama)
  #nparam<-length(stime)
  #   if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
  #     ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
  #     ICOMP<-trvarmat<-NA
  #     errstd<-wald<-pval<-rep(NA,nparam)
  #   } else {
  #     ddd<-diag(sqrt(1/diag(varmat)))
  #     cormat<-(ddd%*%varmat)%*%ddd  
  #     trvarmat<-sum(diag(varmat))
  #     ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
  #     errstd<-sqrt(diag(varmat));  wald<-stime/errstd;
  #     pval<-2*(1-pnorm(abs(wald)))
  #   }
  #   rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  
  ####################################################################
  
  results<-list('estimates'=stime,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBpq,'vmatLouis'=vmatLouis)
  #class(results)<-"cub"
  return(results)
}
