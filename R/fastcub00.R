#' @title Main function for CUB models without covariates
#' @description Function to estimate and validate a CUB model without covariates for given ordinal responses.
#' @aliases fastcub00
#' @usage fastcub00(m,ordinal,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param starting Starting values for the algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @param verbose Logical: should messages about acceleration steps be printed out? (Default is FALSE)
#' @return An object of the class "fastCUB"
#' @seealso \code{\link{probbit}}, \code{\link{probcub00}}
#' @keywords internal 
#' @import stats graphics




fastcub00<-function(m,ordinal,starting=NULL,maxiter,toler,iterc=3,invgen=TRUE,verbose=FALSE){
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
  beta0jj<-log(pai/(1-pai)); gama0jj<-log(csi/(1-csi))
  ##################################################################
  loglik<-loglikcub00(m,freq,beta0jj,gama0jj)
  
  duratait<-c()
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,0) ***************************
  # ********************************************************************
  nniter<-1;   
  
  while(nniter<=maxiter){
    ttit<-proc.time()
    thetaold<-c(beta0jj,gama0jj)
    csi<-1/(1+exp(-gama0jj))
    pai<-1/(1+exp(-beta0jj))
    
    likold<-loglik
    bb<-probbit(m,csi)[ordinal]
    aa<-(1-pai)/(m*pai*bb)
    ttau<-1/(1+aa)
    ft<-ordinal*ttau
    averpo<-(sum(ft))/sum(ttau)
    pai<-(sum(ttau))/n  
    csi<-(m-averpo)/(m-1)   
    if(csi<0.001){
      csi<-0.001;nniter<-maxiter-1;
    }
    
    beta0jj<-log(pai/(1-pai)); gama0jj<-log(csi/(1-csi))
    
    thetanew<-c(beta0jj,gama0jj)
    
    ai<- ordinal - 1 -  (m-1)*(1-csi)

    
    if (nniter>= iterc){
      
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
      
      
      tt1<-proc.time()
      
      inv<-Iinv%*%InfC
      dk<- thetanew-thetaold
      
      thetastar<- thetanew + inv%*%dk
      if (verbose==TRUE) {print("accelerating")}
      
      beta0jj<-thetastar[1]; gama0jj<-thetastar[2]
      
      
    } 
    
    
    loglik<-loglikcub00(m,freq,beta0jj,gama0jj)
    liknew<-loglik
    testll<-abs(liknew-likold) ###### 
    if(testll<=toler) break else {loglik<-liknew}

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
  varianze<-diag(vmatLouis) 
  gama0hat<-log(csi/(1-csi)); beta0hat<-log(pai/(1-pai))
  secsi<-sqrt(varianze[2])*exp(-gama0hat)/((1+exp(-gama0hat))^2)
  sepai<-sqrt(varianze[1])*exp(-beta0hat)/((1+exp(-beta0hat))^2)
  se<-c(sepai,secsi)
  
  
  
  
  ######
  if(csi>0.999) csi<-0.99                                                             
  if(csi<0.001) csi<-0.01        
  if(pai<0.001) pai<-0.01         
  
  ######
  AICCUB00<- -2*loglik+2*(2)
  BICCUB00<- -2*loglik+log(n)*(2)
  
  
  results<-list('estimates'=stime,'se'=se,'time'=durata, 
                'loglik'=loglik,'niter'=nniter,
                'BIC'=BICCUB00,'vmatLouis'=vmatLouis,'duratait'=duratait)
  
  
  return(results)
}