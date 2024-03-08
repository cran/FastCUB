#' @title  S3 method "fitted" for class "fastCUB"
#' @description S3 method fitted for objects of class \code{\link{fastCUB}}. 
#' @aliases fitted.fastCUB
#' @param object An object of class \code{\link{fastCUB}}
#' @param ...  Other arguments
#' @method fitted fastCUB
#' @export 
#' @details Returns the fitted probability distribution for GEM models with no covariates. If only one dichotomous 
#' covariate is included in the model to explain some components, it returns the fitted probability distribution for each profile.
#' @import methods
#' @seealso \code{fastCUB}
#' @rdname fitted.fastCUB
#' @keywords package



fitted.fastCUB<-function(object, ...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }

  theorpr<-profiles(object)
  return(round(theorpr,digits=digits))
}



profiles <- function(object) UseMethod("profiles", object)


profiles.fastCUB<-function(object){
  
  ellipsis<-object$ellipsis
  m<-ellipsis[['m']]
  
  values<-c()
  for (j in 1:m){
    values[j]<-paste("R =",j)
  }
  
  stime<-object$estimates
  
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
 
    if ( is.null(W) & is.null(Y)){
      #  nomi<-rbind("pai","csi");
      pai<-stime[1];csi<-stime[2];
      
      nprof<-1
      theorpr<-matrix(NA,nrow=m,ncol=nprof)
      theorpr[,1]<-probcub00(m,pai,csi)
      profili<-""
      dimnames(theorpr)<-list(values,profili)
      return(theorpr)
      
    } 
    if (is.null(W) & !is.null(Y)) {
      bet<-stime[1:(length(stime)-1)]; csi<-stime[length(stime)]
      #nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),"csi   ")
      #Y<-as.matrix(ellipsis$Y)
      if (NCOL(Y)==1 && length(unique(Y))==2) {
        
        Y<-as.matrix(Y)
        ny<-NCOL(Y)
        yval<-list()
        for (j in 1:ny){
          yval[[j]]<-sort(unique(Y[,j]))
        }
        
        profiles<-expand.grid(yval)
        nprof<-NROW(profiles)
        
        theorpr<-matrix(NA,nrow=m,ncol=nprof)
        
        #paivett<-unique(logis(Y,bet))
        paivett<-c()
        profili<-c()
        for (j in 1:nprof){
          profili[j]<-"("
          paivett[j]<-1/(1+ exp(-bet[1]-sum(as.numeric(profiles[j,])*bet[2:length(bet)])))
          
          theorpr[,j]<-probcub00(m,paivett[j],csi)
          for (k in 1:NCOL(profiles)){
            profili[j]<-paste(profili[j],"Y",k,"=",profiles[j,k],"")
          }
          profili[j]<-paste(profili[j],")")
        }
        dimnames(theorpr)<-list(values,profili)
        return(theorpr)
      } else {
        cat("No fitted method available","\n")
      }
      
    } 
    
    if (!is.null(W) & is.null(Y)){
    #  W<-as.matrix(ellipsis$W)
      
      if (NCOL(W)==1 && length(unique(W))==2){
        pai<-stime[1]; gama<-stime[2:length(stime)];
        # nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"))
        wval<-list()
        
        nw<-NCOL(W)
        W<-as.matrix(W)
        for (j in 1:nw){
          wval[[j]]<-sort(unique(W[,j]))
        }
        profiles<-expand.grid(wval)
        nprof<-NROW(profiles)
        theorpr<-matrix(NA,nrow=m,ncol=nprof)
        csivett<-c()
        #csivett<-unique(logis(W,gama))
        profili<-c()
        for (j in 1:nprof){
          profili[j]<-"("
          csivett[j]<-1/(1+ exp(-gama[1]-sum(as.numeric(profiles[j,])*gama[2:length(gama)])))
          theorpr[,j]<-probcub00(m,pai,csivett[j])
          for (k in 1:NCOL(profiles)){
            profili[j]<-paste(profili[j],"W",k,"=",profiles[j,k])
          }
          profili[j]<-paste(profili[j],")")
        }
        dimnames(theorpr)<-list(values,profili)
        return(theorpr)
        
      }  else {
        cat("No fitted method available","\n")
        
      }
    }
      
    if (!is.null(Y) & !is.null(W)){
     # Y<-as.matrix(ellipsis$Y)
    #  W<-as.matrix(ellipsis$W)
      ny<-NCOL(Y)
      nw<-NCOL(W)
      
      if (ny==1 & nw==1 & length(unique(Y))==2 & length(unique(W))==2 ){
        
        if (all(Y==W)){
          
          bet<-stime[1:(ny+1)];gama<-stime[(ny+2):length(stime)];
          #nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),paste("gamma",0:(length(gama)-1),sep="_"))
          
          listW<-list()
          listY<-list()
          Y<-as.matrix(Y)
          W<-as.matrix(W)
          for (j in 1:ny){
            listY[[j]]<-Y[,j]
          }
          
          for (j in 1:nw){
            listW[[j]]<-W[,j]
          }
          
          eqy<-eqw<-c()
          for (j in 1:length(listY)){
            for (k in 1:length(listW)){
              if (all(listY[[j]]==listW[[k]])){
                eqy<-c(eqy,j)
                eqw<-c(eqw,k)
              }
            }
          }
          comcov<-cbind(eqy,eqw)
          
          YW<-as.matrix(unique(t(cbind(Y,W))))
          #unique restituisce le righe di un array tolte le ripetizioni
          
          paivett<-csivett<-c()
          ywval<-list()
          for (j in 1:NCOL(t(YW))){
            ywval[[j]]<-sort(unique(t(YW)[,j]))
          }
          
          
          profiles<-unique(expand.grid(ywval))
          nprof<-NROW(profiles)
          theorpr<-matrix(NA,nrow=m,ncol=nprof)
          profili<-c()
          # csivett<-unique(logis(W,gama))
          # paivett<-unique(logis(Y,bet))
          
          for (j in 1:nprof){
            vett<-rep(NA,nw)
            if (NROW(comcov)!=0){
              vett[eqw]<-as.numeric(profiles[j,eqy])
              if (nw - NROW(comcov) > 0){
                vett[-eqw]<-as.numeric(profiles[j,(ny+1):NCOL(profiles)])
              }
            } else {
              vett<-as.numeric(profiles[j,(ny+1):NCOL(profiles)])
            }
            profili[j]="("
            for (k in 1:ny){
              profili[j]<-paste(profili[j],"Y",k,"=",profiles[j,k],", ")
            }
            for (k in 1:nw){
              profili[j]<-paste(profili[j],"W",k,"=",vett[k],", ")
            }
            profili[j]<-paste(profili[j],")")
            
            csivett[j]<- 1/(1+ exp(-gama[1]-sum(vett*gama[2:length(gama)])))
            paivett[j]<-1/(1+ exp(-bet[1]-sum(profiles[j,1:ny]*bet[2:length(bet)])))
            theorpr[,j]<-probcub00(m,paivett[j],csivett[j])
          }
          dimnames(theorpr)<-list(values,profili)
          return(theorpr)
          
        } else {
          cat("No fitted method available","\n")
        }
      } else {
        cat("No fitted method available","\n")
      }   
  } 
  
}

