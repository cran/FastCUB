#' @title Main function for CUB models without covariates
#' @description Function to estimate and validate a CUB model without covariates for given ordinal responses.
#' @aliases fastcub00
#' @usage fastcub00(m,ordinal,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @return An object of the class "fastCUB"
#' @seealso \code{\link{probbit}}, \code{\link{probcub00}}
#' @keywords internal 
#' @import stats graphics




fastcub00<-function(m,ordinal,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE){
  tt0<-proc.time()
  serie<-1:m
  freq<-tabulate(ordinal,nbins=m)
  n<-sum(freq)
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  #######################################################
  if (is.null(starting)){
    inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; csi<-inipaicsi[2];
  } else {
    pai<-starting[1]; csi<-starting[2]
  }
  
  ##################################################################
  loglik<-loglikcub00(m,freq,pai,csi)
  
  duratait<-c()
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,0) ***************************
  # ********************************************************************
  nniter<-1;   
  
  while(nniter<=maxiter){
    ttit<-proc.time()
    thetaold<-c(pai,csi)
    
    
    likold<-loglik
    bb<-probbit(m,csi)[ordinal]
    aa<-(1-pai)/(m*pai*bb)
    ttau<-1/(1+aa)
    ft<-ordinal*ttau
    averpo<-(sum(ft))/sum(ttau)
    pai<-(sum(ttau))/n  # updated pai estimate
    csi<-(m-averpo)/(m-1)   # updated csi estimate
    if(csi<0.001){
      csi<-0.001;nniter<-maxiter-1;
    }
    
    
    thetanew<-c(pai,csi)
    
    ai<- ordinal - 1 -  (m-1)*(1-csi)
    
    
    #ai<- (m-1)*(1-csii) - (ordinal-1)
    if (nniter>= iterc){
      
      param<-thetanew
      # M1<-vcScorec(ttau,param,ai,Y=NULL,W=NULL)
      # M2<-vcScore(ttau,param,ai,Y=NULL,W=NULL)
      # 
      # InfC<-Ic(ttau,ordinal,param,ai,Y=NULL,W=NULL)
      bb<-probbit(m,csi)[ordinal]
      aa<-(1-pai)/(m*pai*bb)
      ttau<-1/(1+aa)
      
      dc<-decomp(ttau,ordinal,m,param,ai,Y=NULL,W=NULL)
      
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
      
      
      tt1<-proc.time()
      
      inv<-Iinv%*%InfC
      dk<- thetanew-thetaold
      
      thetastar<- thetanew + inv%*%dk
      
      
      pai<-thetastar[1]; csi<-thetastar[2]
      
      
    } 
    
    # print(c(pai,csi));
    loglik<-loglikcub00(m,freq,pai,csi)
    liknew<-loglik
    testll<-abs(liknew-likold) ###### print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi));
    if(testll<=toler) break else {loglik<-liknew}
    # OPTIONAL printing: print(loglik);
    timeit<-proc.time()-ttit; duratait[nniter]<-timeit[1];
    
    nniter<-nniter+1
  }
  
  if (nniter < iterc){
    param<-thetanew
    bb<-probbit(m,csi)[ordinal]
    aa<-(1-pai)/(m*pai*bb)
    ttau<-1/(1+aa)
    
    dc<-decomp(ttau,ordinal,m,param,ai,Y=NULL,W=NULL)
    
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
  
  
  vmatLouis<-Iinv
  durata<-proc.time()-tt0;durata<-durata[1];
  stime<-c(pai,csi)
  ######
  if(csi>0.999) csi<-0.99                                                             
  if(csi<0.001) csi<-0.01         ### to avoid division by 0 !!!
  if(pai<0.001) pai<-0.01         ### to avoid division by 0 !!!
  ### to ensure identifiability !!!
  ######
  AICCUB00<- -2*loglik+2*(2)
  BICCUB00<- -2*loglik+log(n)*(2)
  
  
  results<-list('estimates'=stime,'time'=durata, 
                'loglik'=loglik,'niter'=nniter,
                'BIC'=BICCUB00,'vmatLouis'=vmatLouis,'duratait'=duratait)
  
  #class(results)<-"cub"
  return(results)
}
