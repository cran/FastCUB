#' @title Main function for CUB models with covariates for the feeling component
#' @description Function to estimate and validate a CUB model for given ordinal responses, with covariates for
#'  explaining the feeling component.
#' @aliases fastcub0q
#' @usage fastcub0q(m,ordinal,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE, verbose=FALSE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @param verbose Logical: should messages about  acceleration steps be printed out? (Default is FALSE)
#' @import stats graphics
#' @return An object of the class "fastCUB"
#' @keywords internal 


#######################
#######################
fastcub0q<-function(m,ordinal,W,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE, verbose=FALSE){
  tt0<-proc.time()
  n<-length(ordinal)
  W<-as.matrix(W)
 
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
  beta0jj<-log(paijj/(1-paijj))
  loglikjj<-loglikcub0q(m,ordinal,W,beta0jj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,q) ***************************
  # ********************************************************************
  nniter<-1
  iterc2<-iterc
  while(nniter<=maxiter){
    
    thetaold<-c(beta0jj,gamajj)
    paijj<-1/(1+exp(-beta0jj))
    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)
    ttau<-1/(1+(1-paijj)/(m*paijj*vettn)) 
    
    ################################# maximize w.r.t. gama  ########
    ordd<-ordinal;covar<-WW;
    gama<-gamajj
    optimgama<-optim(gama,effe01,esterno01=cbind(ttau,ordinal,WW),m=m) 
    ################################################################
    gamajj<-optimgama$par
    paijj<-sum(ttau)/n                   
    beta0jj<-log(paijj/(1-paijj))
    
    ## aggiungo
    thetanew<-c(beta0jj,gamajj)
    param<-thetanew
    loglikjj<-loglikcub0q(m,ordinal,W,beta0jj,gamajj)
    
   
    if (nniter>= iterc2){
      
      
      csii<-logis(W,gamajj)
      ai<- ordinal - 1 -  (m-1)*(1-csii)
      
      vettn<-as.numeric(bitgama(m,ordinal,W,gamajj)   )
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
        Iinv<- try(solve(D+H),silent=TRUE)
      }     
      
      if (class(Iinv)=="try-error"){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if (any(is.infinite(Iinv)) | any(is.nan(Iinv))){
        Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
      }
      
      if ((class(Iinv)!="try-error") && !any(is.na(Iinv))){
        
        if (verbose==TRUE){
          print("accelerating")
        }
   
        inv<-Iinv%*%InfC
        dk<- thetanew-thetaold
        
        thetastar<- thetanew + inv%*%dk
        loglikold2<-loglikcub0q(m,ordinal,W,beta0jj,gamajj)
        beta0jj2<-thetastar[1]; gamajj2<-thetastar[2:(q+2)]
        loglikjj2<-loglikcub0q(m,ordinal,W,beta0jj2,gamajj2)
        
        if ( is.nan(loglikjj2)) {
          Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
         if (verbose==TRUE) {cat("not accelerating at iter",nniter,"\n")}
          
          nniter <- nniter + 1
          next
        } else if (loglikjj2 < loglikold2){
          Iinv<-matrix(NA,nrow=nrow(D),ncol=nrow(D))
          if (verbose ==TRUE) {cat("not accelerating at iter",nniter,"\n")}
          
          nniter <- nniter + 1
          next
        } else {
          beta0jj<-beta0jj2; gamajj<-gamajj2; paijj<-1/(1+exp(-beta0jj))            #%
          loglikold<-loglikold2
          loglikjj<-loglikjj2
        
        }
        
        
      }  else {
        if (verbose == TRUE) {cat("not accelerating at iter",nniter,"\n")}
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
  
  pai<-paijj;  gama<-gamajj;  loglik<-loglikjj; stime<-c(paijj,gamajj)
  if (nniter < iterc) {
    param<-c(beta0jj,gamajj)
    
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
      Iinv<- try(solve(D+H),silent=TRUE)
    }     
  }
  
  if (class(Iinv)=="try-error"){
    Iinv<- matrix(NA,nrow=length(param),ncol=length(param))
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
      varianze<-diag(vmatLouis)
      paihat<-stime[1]
      beta0hat<-log(paihat/(1-paihat))
      sepai<-sqrt(varianze[1])*exp(-beta0hat)/((1+exp(-beta0hat))^2)
      se<-c(sepai,sqrt(varianze[-1]))
      
    }
  }
  durata<-proc.time()-tt0;durata<-durata[1];
  
  
  ####################################################################
  AICCUB0q<- -2*loglik+2*(q+2)
  BICCUB0q<- -2*loglik+log(n)*(q+2)
  ####################################################################
 
  
  
  results<-list('estimates'=stime,'time'=durata,'se'=se,'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=BICCUB0q,'vmatLouis'=vmatLouis)
 
  return(results)
}
