#' @title Main function for CUB models with covariates for the uncertainty component
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining 
#' the feeling component via a logistic transform.
#' @aliases fastcubp0
#' @usage fastcubp0(m,ordinal,Y,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param verbose Logical: should messages about  acceleration steps be printed out? (Default is FALSE)
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @return An object of the class "fastCUB"
#' @import stats graphics
#' @keywords internal 



fastcubp0<-function(m,ordinal,Y,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE){
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
    betjj<- c(bet0,rep(0.1,p))               
    csijj<-inipaicsi[2]
  } else {
    betjj<-starting[1:(p+1)]; csijj<-starting[p+2]
  }
  
  ##############################################################
  gama0jj<-log(csijj/(1-csijj))
  loglikjj<-loglikcubp0(m,ordinal,Y,betjj,gama0jj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,0) ***************************
  # ********************************************************************
  nniter<-1
  iterc2<-iterc
  
  while(nniter<=maxiter){
    
    thetaold<-c(betjj,gama0jj)
    
    loglikold<-loglikjj
    bb<-probbit(m,csijj)
    vettn<-bb[ordinal]     
    aai<- -1+ 1/(logis(Y,betjj)) 
    ttau<-1/(1+aai/(m*vettn))       
    averpo<-sum(ordinal*ttau)/sum(ttau)
    ################################## maximize w.r.t. bet  ########
    bet<-betjj
    covar<-YY
    tauno<-ttau
    
    opmaxbet<-optim(bet,effe10,esterno10=cbind(tauno,covar))
    ################################################################         
    betjj<-opmaxbet$par
   
    csijj<-(m-averpo)/(m-1)       
  
    gama0jj<-log(csijj/(1-csijj))
    loglikjj<-loglikcubp0(m,ordinal,Y,betjj,gama0jj)
    
    
    thetanew<-c(betjj,gama0jj)
    param<-thetanew
    
    if (nniter>= iterc2){
      ai<- ordinal - 1 -  (m-1)*(1-csijj)
      
      bb<-probbit(m,csijj)
      vettn<-bb[ordinal] 
      aai<- -1+1/(logis(Y,betjj))
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
        Iinv<- try(solve(D+H),silent=TRUE)
      } 
      
      if (class(Iinv)=="try-error"){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if (any(is.infinite(Iinv)) | any(is.nan(Iinv))){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if ((class(Iinv)!="try-error") && !any(is.na(Iinv))){
        
        if (verbose==TRUE) {print("accelerating")}
        
        loglikold2<- loglikcubp0(m,ordinal,Y,betjj,gama0jj)
        
        inv<-Iinv%*%InfC
        dk<- thetanew-thetaold
        thetastar<- thetanew + inv%*%dk
        betjj2<-thetastar[1:(p+1)]; gama0jj2<-thetastar[(p+2)]
        loglikjj2<- loglikcubp0(m,ordinal,Y,betjj2,gama0jj2)
        
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
          loglikold<-loglikold2
          loglikjj<-loglikjj2
          betjj<-betjj2; gama0jj<-gama0jj2; csijj<-1/(1+exp(-gama0jj))
        }
        
        

      }  else {
        if (verbose==TRUE) {cat("not accelerating at iter",nniter,"\n")}
       
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
    aai<- -1+ 1/(logis(Y,betjj)) 
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
      Iinv<- try(solve(D+H),silent=TRUE)
    }
    if (class(Iinv)=="try-error"){
      Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
    }
    
    
  }
  
  stime<-c(betjj,csijj)
  bet<-betjj;  csi<-csijj; gama0<-gama0jj; loglik<-loglikjj;
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
      varianze<-diag(varmat)
      p<-NCOL(Y)
      csihat<-stime[-c(1:(p+1))]
      gama0hat<-log(csihat/(1-csihat))
      secsi<-sqrt(varianze[(p+2)])*exp(-gama0hat)/((1+exp(-gama0hat))^2)
      se<-c(sqrt(varianze[-c(p+2)]),secsi)
    }
  }
  
  durata<-proc.time()-tt0;durata<-durata[1];
  
  
  
  
  ####################################################################
  AICCUBp0<- -2*loglik+2*(p+2)
  BICCUBp0<- -2*loglik+log(n)*(p+2)
  ###################################################################
  
  
  
  
  results<-list('estimates'=stime, 'se'=se,  'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBp0,'vmatLouis'=vmatLouis)
  
  return(results)
}
