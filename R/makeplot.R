#' @title Plot facilities for fastCUB objects
#' @description Plot facilities for objects of class "fastCUB". 
#' @aliases makeplot
#' @param object An object of class "fastCUB"
#' @export
#' @details Returns a plot comparing fitted 
#' probabilities and observed relative frequencies for GEM models without covariates. If only one 
#' explanatory dichotomous variable is included in the model for one or all components, 
#' then the function returns a plot comparing the distributions of the responses conditioned to
#' the value of the covariate. 
#' @keywords models device package


makeplot<-function(object){
  
    makeplotCUB(object)
  
}


makeplotCUB<-function(object){
  
  ellipsis<-object$ellipsis
  ordinal<-object$ordinal
  family<-object$family
  m <- ellipsis[['m']]
  
  n<-length(ordinal)
  
  modello<-object$formula
  data<-ellipsis$data
  
  mf<-model.frame(modello,data=data,na.action=na.omit)

  covpai<-model.matrix(modello,data=mf,rhs=1)
  covcsi<-model.matrix(modello,data=mf,rhs=2)

  if (ncol(covpai)==0){
    Y<-NULL
  } else {
    Y<-covpai[,-1]
  }
  if (ncol(covcsi)==0){
    W<-NULL
  } else {
    W<-covcsi[,-1]
  }


  stime<-round(object$estimates,5);

    if ( is.null(W) & is.null(Y) ){
      theorpr<-fitted(object)
      #freq<-matrix(NA,nrow=m,ncol=nprof)
      freq<-tabulate(ordinal,nbins=m)
      stringtitle<-"CUB model";
      thpr<-theorpr[,1]
      dissimi<-dissim(thpr,freq/n)
      #par(mfrow=c(2,1))
      plot(cbind(1:m,1:m),cbind(thpr,(freq/n)),las=1,
           main=paste(stringtitle,  "     (Diss =",round(dissimi,digits=4),")"),
           xlim=c(1,m),ylim=c(0.0,1.1*max(thpr,(freq/n))),
           xlab="Ordinal values of R=1,2,...,m",
           ylab=expression(paste("Observed relative frequencies (dots) and fitted probabilities (circles)")));
      ###
      points(1:m,thpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3); ### ex pch=8,col="red"
      points(1:m,freq/n,pch=16,cex=1.25,lwd=1.5);
      abline(h=0);

      
    } 
      
      if (is.null(W) & !is.null(Y)) {
      

      if (length(unique(Y))==2){
        theorpr<-fitted(object)
        prob0<-theorpr[,1]
        prob1<-theorpr[,2]
        maxpr<-max(prob0,prob1)
        
        plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
             main="CUB distributions, given pai-covariate=0, 1",
             cex=1.2,xlab="Ordinal values of R=1,2,...,m",
             ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
        lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
        abline(h=0);
        
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
        # multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
        
      }
      
    } 
    if (!is.null(W) & is.null(Y) ){
      
      #W<-as.matrix(ellipsis$W)
      
      if (NCOL(W)==1 && length(unique(W))==2){
        theorpr<-fitted(object)
        
        prob0<-theorpr[,1]
        prob1<-theorpr[,2]
        maxpr<-max(prob0,prob1)
        
        plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
             main="CUB distributions, given csi-covariate=0, 1",
             cex=1.2,xlab="Ordinal values of R=1,2,...,m",
             ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
        lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
        abline(h=0);
        
        #  multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
        
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
        

      }
      
    } 
    if(!is.null(W) & !is.null(Y) ){
      

      ny<-NCOL(Y)
      nw<-NCOL(W)
      
      if (ny==1 & nw==1 & length(unique(Y))==2 & length(unique(W))==2 ){
        
        if (all(Y==W)){
          theorpr<-fitted(object)
          #par(mfrow=c(2,1))
          prob0<-theorpr[,1]
          prob1<-theorpr[,2]
          maxpr<-max(prob0,prob1)
          
          plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.lab=0.9,cex.main=0.9,las=1,
               main="CUB distributions, given pai-csi covariate=0, 1",
               cex=1.2,xlab="Ordinal values of R=1,2,...,m",
               ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
          lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
          abline(h=0);
          
    
        }
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
      } 
      
    } 
    
}
