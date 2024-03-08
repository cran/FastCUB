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

bestcub<-
  function (ordinal, m, Y, W, toler = 1e-04, maxiter = 200, iterc = 5, 
            alpha = 0.05, mix = FALSE, tolmix = 100, fmix = NULL, invgen = TRUE) 
  {
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
    vsel <- list()
    matiterfast <- matduratafast <- matrix(nrow = length(setW), 
                                           ncol = length(setY))
    for (i in 1:length(setW)) {
      WW <- as.matrix(W[, setW[[i]]])
      for (j in 1:length(setY)) {
        YY <- as.matrix(Y[, setY[[j]]])
        tt <- (i - 1) * length(setY) + j
        cat("##############################", "\\n")
        cat("Modello", tt, "\\n")
        if (mix == TRUE) {
          cub <- GEM(Formula(ordinal ~ YY | WW | 0), family = "cub", 
                     toler = tolmix)
          iterc <- ifelse(!is.null(fmix), round(fmix * 
                                                  cub$niter), 1)
          starting <- cub$estimates
        }
        else {
          starting = NULL
        }
        fitcov <- fastCUB(Formula(ordinal~YY|WW), m, starting, maxiter, 
                          toler, iterc = iterc, invgen)
        if (mix == TRUE) {
          matiterfast[i, j] <- fitcov$niter + cub$niter
          matduratafast[i, j] <- fitcov$time + cub$time
        }
        else {
          matiterfast[i, j] <- fitcov$niter
          matduratafast[i, j] <- fitcov$time
        }
        vsel[[tt]] <- fitcov
        vsel[[tt]]$nomiW <- colnames(W)[setW[[i]]]
        vsel[[tt]]$nomiY <- colnames(Y)[setY[[j]]]
      }
    }
    bicfast <- as.numeric(lapply(vsel, function(x) x$BIC))
    bic_ord <- sort(bicfast, index = TRUE)
    indici <- bic_ord$ix
    allsig <- 0
    l <- 1
    sigmodel <- c()
    while (allsig == 0) {
      jj <- indici[l]
      lq <- length(vsel[[jj]]$nomiW)
      lp <- length(vsel[[jj]]$nomiY)
      param <- vsel[[jj]]$estimates
      se <- sqrt(diag(vsel[[jj]]$vmatLouis))
      wald <- param/se
      wald2 <- wald[-c(1, lp + 2)]
      if (all(abs(wald2) > qnorm(1 - alpha/2))) {
        allsig <- 1
        sigmodel <- c(sigmodel, jj)
      }
      l <- l + 1
    }
    ll <- l - 1
    bestW <- vsel[[indici[ll]]]$nomiW
    bestY <- vsel[[indici[ll]]]$nomiY
    bestparam <- vsel[[indici[ll]]]$estimates
    bestse <- sqrt(diag(vsel[[indici[ll]]]$vmatLouis))
    bestbic <- bicfast[indici[ll]]
    return(list(vsel = vsel, itercvett = itercvett, bestW = bestW, 
                bestY = bestY, param = bestparam, se = bestse, bic = bestbic, 
                mattime = matduratafast, matiter = matiterfast))
  }




