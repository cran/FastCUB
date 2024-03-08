#' @title Main function for CUB models with covariates for the feeling component
#' @description Function to estimate and validate a CUB model for given ordinal responses, with covariates for
#'  explaining the feeling component.
#' @aliases fastcub0q
#' @usage fastcub0q(m,ordinal,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @import stats graphics
#' @return An object of the class "fastCUB"
#' @keywords internal 


#######################
fastcub0q<-function(m,ordinal,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE){
  tt0<-proc.time()
  n<-length(ordinal)
  W<-as.matrix(W)
  # if (ncol(W)==1){
  #   W<-as.numeric(W)
  # }
  q<-NCOL(W)
  aver<-mean(ordinal)
  WW<-cbind(1,W)                   
  ##############################################################
  freq<-tabulate(ordinal,nbins=m); 
  
  if (is.null(starting)){
    inipaicsi<-inibest(m,freq);  paijj<-inipaicsi[1]; 
    gamajj<-inibestgama(m,ordinal,W)
  } else {
    paijj<-starting[1]; gamajj<-starting[-1]
  }
  
  ##############################################################
  loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,q) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    
    thetaold<-c(paijj,gamajj)
    
    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)
    ttau<-1/(1+(1-paijj)/(m*paijj*vettn)) 
    
    ################################# maximize w.r.t. gama  ########
    ordd<-ordinal;covar<-WW;
    gama<-gamajj
    optimgama<-optim(gama,effe01,esterno01=cbind(ttau,ordinal,WW),m=m) 
    ################################################################
    gamajj<-optimgama$par
    paijj<-sum(ttau)/n                    #updated pai estimate
    
    
    csii<-logis(W,gamajj)
    
    thetanew<-c(paijj,gamajj)
    
    
    ai<- ordinal - 1 -  (m-1)*(1-csii)
    
    #ai<- (m-1)*(1-csii) - (ordinal-1)
    if (nniter>= iterc){
      
      param<-thetanew
      
      loglikold<-loglikcub0q(m,ordinal,W,paijj,gamajj)
      
      # M1<-vcScorec(ttau,param,ai,Y=NULL,W)
      # M2<-vcScore(ttau,param,ai,Y=NULL,W)
      # 
      # InfC<-Ic(ttau,ordinal,param,ai,Y=NULL,W)
      vettn<-bitgama(m,ordinal,W,gamajj)
      ttau<-1/(1+(1-paijj)/(m*paijj*vettn)) 
      dc<-decomp(ttau,ordinal,m,param,ai,Y=NULL,W=W)
      
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
        Iinv<- solve(D+H)
      }     
      
      inv<-Iinv%*%InfC
      dk<- thetanew-thetaold
      
      thetastar<- thetanew + inv%*%dk
      
      paijj<-thetastar[1]; gamajj<-thetastar[2:(q+2)]
      
      
    }
    
    
    loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)## needed for nlm version
    # print(c(nniter,paijj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    
    
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  
  if (nniter < iterc) {
    param<-c(paijj,gamajj)
    
    vettn<-bitgama(m,ordinal,W,gamajj)
    ttau<-1/(1+(1-paijj)/(m*paijj*vettn))
    dc<-decomp(ttau,ordinal,m,param,ai,Y=NULL,W=W)
    
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
      Iinv<- solve(D+H)
    }     
  }
  
  
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
  durata<-proc.time()-tt0;durata<-durata[1];
  
  pai<-paijj;  gama<-gamajj;  loglik<-loglikjj; stime<-c(paijj,gamajj)
  ####################################################################
  AICCUB0q<- -2*loglik+2*(q+2)
  BICCUB0q<- -2*loglik+log(n)*(q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  
  
  results<-list('estimates'=stime,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUB0q,'vmatLouis'=vmatLouis)
  # class(results)<-"cub"
  return(results)
}

