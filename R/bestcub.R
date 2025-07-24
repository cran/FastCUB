#' @title Best-subset variable selection for CUB models via fast EM algorithm
#' @description Perform a best-subset search for CUB models on the basis of the BIC index, by combining all possible covariates'
#'  specification for feeling and for uncertainty parameters 
#' @aliases bestcub
#' @usage bestcub(ordinal,m,Y,W,toler=1e-4,maxiter=200,iterc=5,alpha=0.05,mix=FALSE,
#' tolmix=1e+2,fmix=NULL,invgen=TRUE)
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param Y Matrix of selected covariates for the uncertainty parameter
#' @param W Matrix of selected covariates for the feeling parameter
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param iterc Iteration from which the acceleration strategy starts 
#' @param alpha Significant level for Wald test
#' @param mix Logical: should a first preliminary standard EM be run at toler equal to tolmix? (default is FALSE)
#' @param tolmix Error tolerance for first preliminary EM (if mix=TRUE).
#' @param fmix Fraction of iteration needed for first preliminary EM (if mix=TRUE). Default is null.
#' @param invgen Logical: should the recursive formula for the inverse of the information matrix be considered? (Default is TRUE)
#' @return A list containing the following results: 
#' \item{vsel}{List of all estimated models (with the accelerated EM)}
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

########################################
bestcub<- function (ordinal, m, Y, W, toler = 1e-04, maxiter = 200, iterc = 5, 
                        alpha = 0.05, mix = FALSE, tolmix = 100, fmix = NULL, invgen = TRUE) {
  
  requireNamespace("utils", "combn")
  #subsetting Y
  if (!is.null(Y)){
    sy <- NCOL(Y)
    ss <- 0
    itercvett <- c()
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
  } else {
    setY<-list(); length(setY)<-1
    setY[[1]]<-NA
  }
  
  #subsetting W
  if (!is.null(W)){
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
  }   else {
    setW<-list(); length(setW)<-1
    setW[[1]]<-NA
  } 
  
  
  vsel <- list()
  matiterfast <- matduratafast <- matrix(nrow = length(setW),ncol = length(setY))
  
  
  nmod<-length(setW)*length(setY)
  
  error_vett<-c()
  
  for (i in 1:length(setW)){
    if (!any(is.na(setW[[i]]))){
      WW <- as.matrix(as.matrix(W)[, setW[[i]]])
      colnames(WW)<-colnames(as.matrix(W))[setW[[i]]]
      formulaW<-paste(colnames(WW),collapse="+")
      YW<-WW
      nomiW<-colnames(as.matrix(WW))
      colnames(YW)<-colnames(WW)
    } else {
      YW<-NULL
      nomiW<-NULL
    }
    for (j in 1:length(setY)){
      
      if (!any(is.na(setY[[j]]))){
        YY <- as.matrix(as.matrix(Y)[, setY[[j]]])
        colnames(YY)<-colnames(as.matrix(Y))[setY[[j]]]
        formulaY<-paste(colnames(YY),collapse="+")
        nomiY<-colnames(as.matrix(YY))
        
        if (!is.null(YW)){
          YW<-cbind(YY,YW)
        } else {
          YW<-as.matrix(YY)
        }
      } else {
        nomiY<-NULL
      }
      
      YW_ok<-as.matrix(YW[,match(unique(colnames(YW)),colnames(YW))])
      colnames(YW_ok)<-colnames(YW)[match(unique(colnames(YW)),colnames(YW))]
      dd<-data.frame(ordinal,YW_ok)
      tt<- length(setY)*(i-1) + j      
      
      vsel[[tt]]<-list()
      cat("##############################", "\n")
      cat("Estimating Model n.", tt, "\n")
      
      if (!is.null(nomiW)){
        formulaW<-paste(nomiW,collapse="+")
        if (!is.null(nomiY)){
          formulaY<-paste(nomiY,collapse="+")
        } else {
          formulaY<-0
        }
      } else {
        formulaW<-0
      }
      
      if (mix == TRUE) {
        effe_ini<-paste(formulaY,formulaW,sep="|")
        effe_ini<-paste(effe_ini,0,sep="|")
        cub <- GEM(as.Formula(paste("ordinal",effe_ini,sep="~")), family = "cub",toler = tolmix,data=dd,m=m)
        iterc <- ifelse(!is.null(fmix), round(fmix*cub$niter), 1)
        starting <- cub$estimates
      }
      else {
        starting = NULL
      }
      
      fitcov <- try(fastCUB(as.Formula(paste("ordinal",paste(formulaY,formulaW,sep="|"),sep="~")), data=dd,m, starting, maxiter, 
                            toler, iterc = iterc, invgen=invgen),silent=TRUE)
      
      if (any(class(fitcov)=="try-error")){
        matiterfast[i, j] <- NA
        matduratafast[i, j] <- NA
        print("error")
        error_vett<-c(error_vett,tt)
        #print(error_vett)
        
      } else {
        if (mix == TRUE) {
          matiterfast[i, j] <- fitcov$niter + cub$niter
          matduratafast[i, j] <- fitcov$time + cub$time
        }
        else {
          matiterfast[i, j] <- fitcov$niter
          matduratafast[i, j] <- fitcov$time
        }
      }
      #names(vsel[[tt]])<-tt
      vsel[[tt]][[1]] <- fitcov

      vsel[[tt]][[2]] <- nomiW 
      vsel[[tt]][[3]] <- nomiY
      names(vsel[[tt]])<-c('model','nomiW','nomiY')
      # if (class(vsel[[tt]]$model)=="try-error"){
      #   
      #   # vsel[[tt]]$model<-list()
      #   # vsel[[tt]]$model$BIC<-Inf
      #   # vsel[[tt]]$model$estimates<-NA
      #   # vsel[[tt]]$model$varmat<-NA
      # }
      
    } #fine ciclo j
  }# fine ciclo i
  
  
  if (!is.null(error_vett)){
    vsel2<-list()
    vsel2<-vsel[-error_vett]
    names(vsel2)<-names(vsel)[-error_vett]
    vsel<-vsel2
    
  }
  bicvett<-lapply(vsel, function(x) x[[1]]$BIC)
  bicfast <- as.numeric(bicvett)
  
  bic_ord <- sort(bicfast, index = TRUE)
  indici <- bic_ord$ix
  allsig <- 0; nosig<-0
  l <- 1
  sigmodel <- c()
  while ((allsig == 0) & (l <= length(indici))) {
    jj <- indici[l]
    if (!is.null(vsel[[jj]]$nomiW)){
      lq <- length(vsel[[jj]]$nomiW)
    } else {
      lq<-0
    }
    if (!is.null(vsel[[jj]]$nomiY)){
      lp <- length(vsel[[jj]]$nomiY)
    } else {
      lp<-0
    }
    
    param <- vsel[[jj]]$model$estimates
    se <- sqrt(diag(vsel[[jj]]$model$varmat))
    wald <- param/se
    wald2 <- wald[-c(1, lp + 2)]
    if (all(abs(wald2) > qnorm(1 - alpha/2))) {
      allsig <- 1
      sigmodel <- c(sigmodel, jj)
    }
    l<-l+1
  }
  if (allsig!=0){
    ll <- l - 1
  } else {
    ll<- 1
    print("Warning: No significant model at selected alpha level")
    
  }
  bestW <- vsel[[indici[ll]]]$nomiW
  bestY <- vsel[[indici[ll]]]$nomiY
  bestparam <- vsel[[indici[ll]]]$model$estimates
  bestse <- sqrt(diag(vsel[[indici[ll]]]$model$varmat))
  bestbic <- bicfast[indici[ll]]
  
  best_model<-vsel[[indici[ll]]]
  
  return(list(vsel = vsel, bestW = bestW, 
              bestY = bestY, param = bestparam, se = bestse, bic = bestbic, 
              mattime = matduratafast, matiter = matiterfast,'bicfast'=bicfast,'best_model'=best_model))
  
}

