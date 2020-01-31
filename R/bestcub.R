#' @title Best-subset variable selection for CUB models via fast EM algorithm
#' @description Perform a best-subset search for CUB models on the basis of the BIC index, by combining all possible covariates'
#'  specification for feeling and for uncertainty parameters
#' @aliases bestcub
#' @usage bestcub(ordinal,m,Y,W,toler=1e-4,maxiter=200,iterc=5,alpha=0.05,mix=FALSE,
#' tolmix=1e+2,fmix=NULL,invgen=TRUE,verbose=FALSE)
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param Y Matrix of selected covariates for the uncertainty parameter
#' @param W Matrix of selected covariates for the feeling parameter
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm
#' @param toler Fixed convergence tolerance for final estimates
#' @param iterc Iteration from which the acceleration strategy starts
#' @param alpha Significant level for Wald test
#' @param mix Logical: should a first preliminary standard EM be run at toler equal to tolmix? (default is FALSE)
#' @param tolmix Convergence tolerance for first preliminary EM (if mix=TRUE).
#' @param fmix Fraction of iteration needed for first preliminary EM (if mix=TRUE). Default is null.
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @param verbose Logical: should messages about  acceleration steps be printed out? (Default is FALSE)
#' @return A list containing the following results:
#' \item{vsel}{List of all estimated models (with the accelerated EM) as FastCUB objects}
#' \item{bestmodel}{FastCUB object corresponding to the best CUB model (according to BIC), if not null}
#' \item{bestW}{Names of covariates for feeling in the best model with all significant effect}
#' \item{bestY}{Names of covariates for feeling in the best model with all significant effect}
#' \item{param}{ML estimates of the best model}
#' \item{se}{Estimated standard errors for the best model}
#' \item{bic}{BIC index of the best (significant) model}
#' \item{mattime}{Matrix of computational time for each of the estimated model}
#' \item{matiter}{Matrix of number of iterations occurred for each of the estimated model}
#' @importFrom CUB GEM
#' @importFrom utils combn
#' @export bestcub
#' @seealso \code{\link{fastCUB}}
#' @keywords stats
#' @examples
#' \donttest{
#' library(FastCUB)
#' data(univer)
#' ordinal<-univer$global
#' m<-7
#' Y<-univer[,c(2,3,4)]
#' W<-univer[,8:11]
#' ## Search for the best CUB model with covariates only for feeling
#' best0q<-bestcub(ordinal,m,Y=NULL,W,toler=1e-4,maxiter=100,iterc=5,alpha=0.05,invgen=TRUE)
#' ## Search for the best CUB model with covariates only for uncertainty
#' bestp0<-bestcub(ordinal,m,Y,W=NULL,toler=1e-4,maxiter=100,iterc=5,alpha=0.05,invgen=TRUE)
#' ## Search for the best CUB model with covariates for both parameters
#' bestpq<-bestcub(ordinal,m,Y,W,toler=1e-4,maxiter=100,iterc=5,alpha=0.05,invgen=TRUE,
#'     mix=TRUE,tolmix=1e+3,fmix=1)
#' final<-bestpq$bestmodel; summary(final)
#' }
#'



########################################


bestcub<- function (ordinal, m, Y, W, toler = 1e-04, maxiter = 200, iterc = 5,
                    alpha = 0.05, mix = FALSE, tolmix = 100, fmix = NULL, invgen = TRUE,verbose=FALSE)
{

  if (!is.null(Y)){
    Y<-as.matrix(Y)
    p<-NCOL(Y)
    if (p>0){
      if (is.null(colnames(Y))){
        colnames(Y)<-rep(NA,p)
        for (h in 1:p){
          colnames(Y)[h]<-paste("Y",h,sep="_")
        }
      }
    }
  }
  if (!is.null(W)){
    W<-as.matrix(W)
    q<-NCOL(W)
    if (q>0){
      if (is.null(colnames(W))){
        colnames(W)<-rep(NA,q)
        for (h in 1:q){
          colnames(W)[h]<-paste("W",h,sep="_")
        }
      }
    }
  }




  sy <- NCOL(Y)
  ss <- 0
  itercvett <- c()
  requireNamespace("utils", "combn")
  setY <- list()
  for (k in 1:sy) {
    indy <- combn(sy, k)
    siy <- NCOL(indy)
    for (l in 1:siy) {
      indyy <- indy[, l]
      ss <- ss + 1
      setY[[ss]] <- indyy
    }
  }
  sw <- NCOL(W)
  ss <- 0
  setW <- list()
  for (k in 1:sw) {
    indw <- combn(sw, k)
    siw <- NCOL(indw)
    for (l in 1:siw) {
      indww <- indw[, l]
      ss <- ss + 1
      setW[[ss]] <- indww
    }
  }
  vsel <- modelli<-list()
  isna<-c()
  matiterfast <- matduratafast <- matrix(nrow = length(setW), ncol = length(setY))
  for (i in 1:length(setW)){
    WW <- W[, setW[[i]]]
    if (NCOL(WW)==1 &  !is.null(WW)){
      WW<-as.matrix(WW)
      colnames(WW)<-colnames(W)[setW[[i]]]
    }
    for (j in 1:length(setY)){
      YY <- Y[, setY[[j]]]
       if (NCOL(YY)==1 & !is.null(YY)){
         YY<-as.matrix(YY)
         colnames(YY)<-colnames(Y)[setY[[j]]]
       }
      tt <- (i - 1)*length(setY) + j
      if (verbose==TRUE){
        cat("##############################", "\n")
        cat("Model N.", tt, "\n")
        print("Covariate csi:")
        print(colnames(WW))
        print("Covariate pai:")
        print(colnames(YY))
      }
      
      
      vsel[[tt]]<-list()
      if (!is.null(YY) & !is.null(WW)){
        effe<-Formula(ordinal ~ YY|WW)
        effecub<-Formula(ordinal ~ YY|WW|0)
      } else {
        if (is.null(YY)){
          effe<-Formula(ordinal ~ 0|WW)
          effecub<-Formula(ordinal ~ 0|WW|0)
        } else {
          effe<-Formula(ordinal ~ YY|0)
          effecub<-Formula(ordinal ~ YY|0|0)
        }
      }



      if (mix == TRUE) {
        requireNamespace("CUB", "GEM")

        cub <- try(GEM(effecub, family = "cub", toler = tolmix),silent=TRUE)
        if (all(class(cub)!="try-error")){
          iterc <- ifelse(!is.null(fmix), round(fmix*cub$niter), 1)
          starting <- cub$estimates
          nitercub<-cub$niter; timecub<-cub$time
        } else {
          starting = NULL
          iterc<-10
          nitercub<-0; timecub<-0
        }
      }  else {
        starting = NULL
      }

    

      dataa <-  environment(effe) 
      fitcov <- fastCUB(effe,dataa, m, starting, maxiter, toler, iterc, invgen,verbose)
     

      if (mix == TRUE) {
        matiterfast[i, j] <- fitcov$niter + nitercub
        matduratafast[i, j] <- fitcov$time + timecub
      }
      else {
        matiterfast[i, j] <- fitcov$niter
        matduratafast[i, j] <- fitcov$time
      }
      
      if (any(class(fitcov)=="try-error")){
        isna<-c(isna,tt)
      }
      
      vsel[[tt]]$modello <- fitcov
      modelli[[tt]]<-fitcov
      vsel[[tt]]$nomiW <- colnames(W)[setW[[i]]]
      vsel[[tt]]$nomiY <- colnames(Y)[setY[[j]]]
    }
  }



  bicfast <- as.numeric(lapply(modelli, function(x) x$BIC))
  
  bic_ord <- sort(bicfast, index = TRUE)
 
  indici <- bic_ord$ix
  allsig <- 0
  l <- 1
  sigmodel <- c()
  while (allsig == 0) {
    jj <- indici[l]
    lq <- length(vsel[[jj]]$nomiW);
    lp <- length(vsel[[jj]]$nomiY);
    param <- vsel[[jj]]$modello$estimates
    ser <-   vsel[[jj]]$modello$se    
    
    wald <- param/ser;
    
    wald2 <- wald[-c(1, lp + 2)]
    
    if (!any(is.na(wald2))){
      
      
      if (all(abs(wald2) > qnorm(1 - alpha/2))) {
        allsig <- 1
        sigmodel <- c(sigmodel, jj); 
        
      }
    }
    l<-l+1; 
  }
  ll <- l-1

  if (is.na(sigmodel)){
    message("No significant model at given alpha level")
    return(list(vsel = vsel, mattime = matduratafast, matiter = matiterfast))

  } else {

    final<-  modelli[[indici[ll]]]
   
    bestW <- vsel[[indici[ll]]]$nomiW
    bestY <- vsel[[indici[ll]]]$nomiY
    bestparam <- final$estimates
    bestse <- final$se
    bestbic <- bicfast[indici[ll]]
    

     nomipai<-ifelse(!is.null(bestY),"intercept","pai")
     nomicsi<-ifelse(!is.null(bestW),"intercept","csi")
    
    listanomi<-c(nomipai,bestY,nomicsi,bestW)
   
   if (verbose==TRUE){
     print("Best model:")
     summary(final,pnames=listanomi)
   }
    
   final$parnames<-listanomi
     
    return(list(vsel = vsel, bestmodel=final, bestW = bestW,
              bestY = bestY, param = bestparam, se = bestse, bic = bestbic,
              mattime = matduratafast, matiter = matiterfast))
  }

}


##########



