#' @title Main function for CUB models with covariates for the uncertainty component
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining 
#' the feeling component via a logistic transform.
#' @aliases fastcubp0
#' @usage fastcubp0(m,ordinal,Y,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @return An object of the class "fastCUB"
#' @import stats graphics
#' @keywords internal 


fastcubp0<-function(m,ordinal,Y,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE){
  tt0<-proc.time()
  n<-length(ordinal)
  Y<-as.matrix(Y)
  
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  
  p<-NCOL(Y)
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  YY<-cbind(1,Y)
  ##################################################################
  serie<-1:m; freq<-tabulate(ordinal,nbins=m);
  
  if (is.null(starting)){
    inipaicsi<-inibest(m,freq)
    pai<-inipaicsi[1]; bet0<-log(pai/(1-pai));  
    betjj<- c(bet0,rep(0.1,p))                #betjj<-rep(0.1,p+1);
    csijj<-inipaicsi[2]
  } else {
    betjj<-starting[1:(p+1)]; csijj<-starting[p+2]
  }
  
  ##############################################################
  loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,0) ***************************
  # ********************************************************************
  nniter<-1
  
  
  while(nniter<=maxiter){
    
    thetaold<-c(betjj,csijj)
    
    loglikold<-loglikjj
    bb<-probbit(m,csijj)
    vettn<-bb[ordinal]      # probbit for all ordinal (r_i,i=1,2,...,n)
    aai<- -1+ 1/(logis(Y,betjj)) #exp(-(YY%*%betjj));
    ttau<-1/(1+aai/(m*vettn))       # tau is a reserved word in R
    averpo<-sum(ordinal*ttau)/sum(ttau)
    ################################## maximize w.r.t. bet  ########
    bet<-betjj
    covar<-YY
    tauno<-ttau
    #nlmaxbet<-nlm(effe10,betjj,esterno10);   
    opmaxbet<-optim(bet,effe10,esterno10=cbind(tauno,covar))
    ################################################################         
    betjj<-opmaxbet$par
    # betjj<-nlmaxbet$estimate;        #updated bet estimates
    csijj<-(m-averpo)/(m-1)       #updated csi estimate
    #loglikjj<- -opmaxbet$value
    
    
    thetanew<-c(betjj,csijj)
    
    ai<- ordinal - 1 -  (m-1)*(1-csijj)
    
    #ai<- (m-1)*(1-csii) - (ordinal-1)
    if (nniter>= iterc){
      
      param<-thetanew
      bb<-probbit(m,csijj)
      vettn<-bb[ordinal]      # probbit for all ordinal (r_i,i=1,2,...,n)
      aai<- -1+ 1/(logis(Y,betjj)) #exp(-(YY%*%betjj));
      ttau<-1/(1+aai/(m*vettn))
      ai<- ordinal - 1 -  (m-1)*(1-csijj)
      loglikold<-loglikcubp0(m,ordinal,Y,betjj,csijj)
      
      dc<-decomp(ttau,ordinal,m,param,ai,Y=Y,W=NULL)
      
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
      
      betjj<-thetastar[1:(p+1)]; csijj<-thetastar[p+2]
      
    } 
    loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
    
    #print(c(nniter,betjj,csijj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  
  if (nniter < iterc ){
    param<-thetanew
    aai<- -1+ 1/(logis(Y,betjj)) #exp(-(YY%*%betjj));
    ttau<-1/(1+aai/(m*vettn)) 
    
    dc<-decomp(ttau,ordinal,m,param,ai,Y=Y,W=NULL)
    
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
  
  #vmatLouis<-Iinv 
  
  bet<-betjj;  csi<-csijj;  loglik<-loglikjj;
  ####################################################################
  AICCUBp0<- -2*loglik+2*(p+2)
  BICCUBp0<- -2*loglik+log(n)*(p+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  
  stime<-c(bet,csi)
  
  
  results<-list('estimates'=stime,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBp0,'vmatLouis'=vmatLouis)
  #class(results)<-"cub"
  return(results)
}
